#include <iostream>

#include "admm.hpp"
#include "problem_data/rand_prob_tinympc_params.hpp"
#include "problem_data/rand_prob_tinympc_xbar.hpp"
#include "problem_data/rand_prob_tinympc_delta.hpp"

#include "Arduino.h"

// #include <vector>
// #define DT 1 / 200

extern "C"
{
  union BtoF
  {
    byte b[8];     // 8바이트: 두 개의 4바이트 float 값
    float fval[2]; // 두 개의 float 값을 저장할 배열
  } r;

  const int buffer_size = 8; // float 두 개에 필요한 8바이트
  byte buf[buffer_size];

  float beta, yaw_rate; // 두 개의 float 값을 저장할 변수
  float sumVal;
  int cnt; // step

  // float yaw_moment;
  Eigen::Map<const Matrix<tinytype, NTOTAL, NSTATES, Eigen::RowMajor>> Xref_total(Xref_data);

  void readFromSimulink()
  {
    Serial.readBytes(buf, buffer_size); // 시리얼로부터 8바이트 읽기
    // 버퍼 데이터를 union의 b 배열에 복사
    for (int i = 0; i < buffer_size; i++)
    {
      r.b[i] = buf[i];
    }
    // union에서 두 개의 float 값을 추출
    beta = r.fval[0];
    yaw_rate = r.fval[1];
  }

  void writeToSimulink(float number1, float number2)
  {
    byte buf[8]; // 두 개의 float 값을 담을 바이트 배열

    memcpy(&buf[0], &number1, sizeof(float)); // 첫 번째 float 값을 배열에 복사
    memcpy(&buf[4], &number2, sizeof(float)); // 두 번째 float 값을 배열에 복사

    Serial.write(buf, 8); // 두 개의 float 값을 8바이트로 전송
  }

  void clearSerialBuffer()
  {
    while (Serial.available() > 0)
    {
      Serial.read(); // 버퍼에 남아 있는 데이터를 읽어 버퍼를 비움
    }
  }
  static TinyCache cache;
  static TinyWorkspace work;
  static TinySettings settings;
  static TinySolver solver{&settings, &cache, &work};

  void setup()
  {
    Serial.begin(115200);
    
    Serial.println("Serial initialized");

    // Copy data from problem_data/rand_prob_tinympc_*.hpp
    cache.rho = rho_value;
    cache.Kinf = Eigen::Map<const Matrix<tinytype, NINPUTS, NSTATES, Eigen::RowMajor>>(Kinf_data);
    cache.Pinf = Eigen::Map<const Matrix<tinytype, NSTATES, NSTATES, Eigen::RowMajor>>(Pinf_data);
    cache.Quu_inv = Eigen::Map<const Matrix<tinytype, NINPUTS, NINPUTS, Eigen::RowMajor>>(Quu_inv_data);
    cache.AmBKt = Eigen::Map<const Matrix<tinytype, NSTATES, NSTATES, Eigen::RowMajor>>(AmBKt_data);
    cache.coeff_d2p = Eigen::Map<const Matrix<tinytype, NSTATES, NINPUTS>>(coeff_d2p_data);

    work.Adyn = Eigen::Map<const Matrix<tinytype, NSTATES, NSTATES, Eigen::RowMajor>>(Adyn_data);
    // work.Bdyn = Eigen::Map<const Matrix<tinytype, NSTATES, NINPUTS, Eigen::RowMajor>>(Bdyn_data);
    work.Bdyn = Eigen::Map<const Matrix<tinytype, NSTATES, NINPUTS>>(Bdyn_data);

    work.Q = Eigen::Map<const tiny_VectorNx>(Q_data);
    work.Qf = Eigen::Map<const tiny_VectorNx>(Qf_data);
    work.R = Eigen::Map<const tiny_VectorNu>(R_data);
    work.u_min = tiny_MatrixNuNhm1::Constant(-800.0);
    work.u_max = tiny_MatrixNuNhm1::Constant(800.0);
    // work.x_min = tiny_MatrixNxNh::Constant(-0.5);
    // work.x_max = tiny_MatrixNxNh::Constant(0.5);

    tiny_VectorNx x_min_one_time_step(-0.1745, -0.5);
    tiny_VectorNx x_max_one_time_step(0.1745, 0.5);
    work.x_min = x_min_one_time_step.replicate(1, NHORIZON);
    work.x_max = x_max_one_time_step.replicate(1, NHORIZON);

    // for (int i=0; i<NHORIZON; i++) {
    // work.x_min[i] = tiny_VectorNc::Constant(-99999); // Currently unused
    // work.x_max[i] = tiny_VectorNc::Zero();
    // work.x_min[i] = tiny_MatrixNxNh::Constant(-10000);
    // work.x_max[i] = tiny_MatrixNxNh::Constant(10000);
    // work.A_constraints[i] = tiny_MatrixNcNx::Zero();
    // }
    work.Xref = tiny_MatrixNxNh::Zero();
    work.Uref = tiny_MatrixNuNhm1::Zero();
    // std::cout << "Adyn: " << work.Adyn << std::endl;

    work.x = tiny_MatrixNxNh::Zero();
    work.q = tiny_MatrixNxNh::Zero();
    work.p = tiny_MatrixNxNh::Zero();
    work.v = tiny_MatrixNxNh::Zero();
    work.vnew = tiny_MatrixNxNh::Zero();
    work.g = tiny_MatrixNxNh::Zero();

    work.u = tiny_MatrixNuNhm1::Zero();
    work.r = tiny_MatrixNuNhm1::Zero();
    work.d = tiny_MatrixNuNhm1::Zero();
    work.z = tiny_MatrixNuNhm1::Zero();
    work.znew = tiny_MatrixNuNhm1::Zero();
    work.y = tiny_MatrixNuNhm1::Zero();

    work.primal_residual_state = 0;
    work.primal_residual_input = 0;
    work.dual_residual_state = 0;
    work.dual_residual_input = 0;
    work.status = 0;
    work.iter = 0;

    settings.abs_pri_tol = 0.01;
    settings.abs_dua_tol = 0.01;
    settings.max_iter = 1000;
    settings.check_termination = 1;
    settings.en_input_bound = 1;
    settings.en_state_bound = 1;
    std::cout << "Setup complete" << std::endl;

    float state_E[2] = {0.00847469954733655, 0.0937660177773983};

    // for (int i = 0; i < NTOTAL - NHORIZON; i++)
    // {
    //   delay(10); // 없으면 터미널 이상하게 출력되 는 듯

    //   // Xref_total = Eigen::Map<Matrix<Eigen::RowMajor, NTOTAL, NSTATES, tinytype>>(Xref_data).transpose();
    //   // work.Xref.block<NSTATES, NHORIZON>(0, 0)= Xref_total.transpose().block<NSTATES, NHORIZON>(0, i);
    //   work.Xref.block<NSTATES, NHORIZON>(0, 0) = Xref_total.transpose().block<NSTATES, NHORIZON>(0, i); // step에 맞는 reference 넣어주기

    //   // Serial.print(availableMemory());

    //   unsigned long start = micros();
    //   tiny_solve(&solver);
    //   unsigned long end = micros();
    //   unsigned long elapsed_time = end - start;      // 경과 시간 계산
    //   float elapsed_time_ms = elapsed_time / 1000.0; // 마이크로초를 밀리초로 변환
    //   // printf("%10.4f %10d %10d\n", (work.x.col(0) - work.Xref.col(0)).squaredNorm(), work.iter, end - start);
    //   // // Serial.println(end-start);
    //   // // std::cout << "TIMING: " << end-start << std::endl;
    //   // // Serial.println((work.x.col(0) - work.Xref.col(0)).squaredNorm());

    //   // float squared_difference[2] = {std::pow(work.x.col(0)(0) - work.Xref.col(0)(0), 2), std::pow(work.x.col(0)(1) - work.Xref.col(0)(1), 2)};
    //   // std::cout << "beta_error_squared: " << squared_difference[0] << std::endl;
    //   // std::cout << "yaw_rate_error_squared: " << squared_difference[1] << std::endl;
    //   // std::cout << "u(0): " << work.u.col(0) << std::endl;

    //   float beta = work.x.col(0)(0);
    //   float yaw = work.x.col(0)(1);

    //   // std::cout << "time: " << i * 0.0005 << std::endl;
    //   // std::cout << "beta: " << beta << std::endl;
    //   // std::cout << "yaw_rate: " << yaw << std::endl;
    //   // std::cout << "u: " << work.u.col(0) << std::endl;

    //   Serial.print(i * 0.01, 6); // time
    //   Serial.print(" ");
    //   Serial.print(beta, 6); // 소수점 이하 6자리까지 출력
    //   Serial.print(" ");
    //   Serial.print(yaw, 6);
    //   Serial.print(" ");
    //   Serial.print(work.u(0, 0), 6);
    //   Serial.print(" ");
    //   Serial.print(elapsed_time_ms);
    //   Serial.print(" ");
    //   Serial.print(work.iter);
    //   Serial.println(); // 새로운 줄로 이동

    //   // writeToSimulink(beta, yaw); // 결과와 경과 시간을 Simulink로 전송
    //   // clearSerialBuffer();                          // clearbuffer 반드시 있어야함.
    //   // delay(100);                                   // delay 반드시 있어야함.

    //   tiny_VectorNx random_num;

    //   for (int j = 0; j < NSTATES; ++j)
    //   {
    //     random_num(j) = state_E[j] * delta[i];
    //   }
    //   tiny_VectorNx next_x = work.Adyn * work.x.col(0) + work.Bdyn * work.u.col(0);
    //   work.x.col(0) = next_x + random_num;
    //   // work.x.col(0) = next_x + E * delta;
    // }
  }
}

void loop()
{
  // delay(1);

  if (Serial.available() >= buffer_size)
  {
    readFromSimulink(); // 데이터 수신

    cnt += 1;
    if (cnt >= (NTOTAL - NHORIZON))
      exit(0);

    // 2. step에 맞는 reference 넣어주기
    work.Xref.block<NSTATES, NHORIZON>(0, 0) = Xref_total.transpose().block<NSTATES, NHORIZON>(0, cnt);
    work.x.col(0) << beta, yaw_rate; // state 넣어주기
    // 3. Solve MPC problem
    //   tiny_solve(&solver);
    //   tiny_VectorNx next_x = work.Adyn * work.x.col(0) + work.Bdyn * work.u.col(0);
    // work.x.col(0) = next_x;

    unsigned long start = micros();
    tiny_solve(&solver); // 연산
    unsigned long end = micros();

    unsigned long elapsed_time = end - start;      // 경과 시간 계산
    float elapsed_time_ms = elapsed_time / 1000.0; // 마이크로초를 밀리초로 변환

    // Serial.println(end-start);

    // std::cout << "TIMING: " << end-start << std::endl;
    //  printf("%10.4f %10d %10d\n", (work.x.col(0) - work.Xref.col(0)).squaredNorm(), work.iter, end - start);

    float yaw_moment = -work.u.col(0)(0);

    writeToSimulink(yaw_moment, elapsed_time_ms); // 결과와 경과 시간을 Simulink로 전송
    clearSerialBuffer();                          // clearbuffer 반드시 있어야함.
    delay(200);                                   // delay 반드시 있어야함.
  }
}

// jh_branch 
// jh_branch 