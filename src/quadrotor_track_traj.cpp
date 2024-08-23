#include <iostream>

#include "admm.hpp"
#include "problem_data/rand_prob_tinympc_params.hpp"
#include "problem_data/rand_prob_tinympc_xbar.hpp"

#include "Arduino.h"

#define DT 1 / 100

extern "C"
{
union BtoF {
  byte b[8];        // 8바이트: 두 개의 4바이트 float 값
  float fval[2];    // 두 개의 float 값을 저장할 배열
} u;

const int buffer_size = 8;  // float 두 개에 필요한 8바이트
byte buf[buffer_size];

float beta, yaw_rate;  // 두 개의 float 값을 저장할 변수
float sumVal;
// float yaw_moment;
Eigen::Map<const Matrix<tinytype, NTOTAL, NSTATES, Eigen::RowMajor>> Xref_total(Xref_data);


void readFromSimulink() {
  Serial.readBytes(buf, buffer_size);  // 시리얼로부터 8바이트 읽기  
  // 버퍼 데이터를 union의 b 배열에 복사
  for (int i = 0; i < buffer_size; i++) {
    u.b[i] = buf[i];
  }
  // union에서 두 개의 float 값을 추출
  beta = u.fval[0];
  yaw_rate = u.fval[1];
}

void writeToSimulink(float number) {
  byte *b = (byte *)&number;
  Serial.write(b, 4);  // 결과 float 값을 4바이트로 전송
  // Serial.write("\r\n");  // CR + LF로 데이터의 끝을 표시
}

void clearSerialBuffer() {
  while (Serial.available() > 0) {
    Serial.read();  // 버퍼에 남아 있는 데이터를 읽어 버퍼를 비움
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
        // cache.coeff_d2p = Eigen::Map<const Matrix<tinytype, NSTATES, NINPUTS, Eigen::RowMajor>>(coeff_d2p_data);

        work.Adyn = Eigen::Map<const Matrix<tinytype, NSTATES, NSTATES, Eigen::RowMajor>>(Adyn_data);
        // work.Bdyn = Eigen::Map<const Matrix<tinytype, NSTATES, NINPUTS, Eigen::RowMajor>>(Bdyn_data);   
        work.Bdyn = Eigen::Map<const Matrix<tinytype, NSTATES, NINPUTS >>(Bdyn_data);        
     
        work.Q = Eigen::Map<const tiny_VectorNx>(Q_data);
        work.Qf = Eigen::Map<const tiny_VectorNx>(Qf_data);
        work.R = Eigen::Map<const tiny_VectorNu>(R_data);
        work.u_min = tiny_MatrixNuNhm1::Constant(-30000.0);
        work.u_max = tiny_MatrixNuNhm1::Constant(30000.0);
        work.x_min = tiny_MatrixNxNh::Constant(-10000);
        work.x_max = tiny_MatrixNxNh::Constant(10000);
        // for (int i=0; i<NHORIZON; i++) {
        // work.x_min[i] = tiny_VectorNc::Constant(-99999); // Currently unused
        // work.x_max[i] = tiny_VectorNc::Zero();
        // work.x_min[i] = tiny_MatrixNxNh::Constant(-10000);
        // work.x_max[i] = tiny_MatrixNxNh::Constant(10000);
        // work.A_constraints[i] = tiny_MatrixNcNx::Zero();
        // }
        work.Xref = tiny_MatrixNxNh::Zero();
        work.Uref = tiny_MatrixNuNhm1::Zero();
        std::cout << "Adyn: " << work.Adyn << std::endl;
       
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

        settings.abs_pri_tol = 0.0001;
        settings.abs_dua_tol = 0.0001;
        settings.max_iter = 4000;
        settings.check_termination = 1;
        settings.en_input_bound = 1;
        settings.en_state_bound = 0;
        std::cout << "Setup complete" << std::endl;
            

    }
 
}

void loop()
{
  if (Serial.available() >= buffer_size) {  
    readFromSimulink(); // 데이터 수신 

    int cnt; //step 
    cnt +=1;
    if(cnt >= (NTOTAL - NHORIZON)) exit(0); 
    work.Xref.block<NSTATES, NHORIZON>(0, 0) = Xref_total.transpose().block<NSTATES, NHORIZON>(0, cnt); //step에 맞는 reference 넣어주기 
    work.x.col(0) << beta, yaw_rate; //state 넣어주기 

    // unsigned long start = micros();
    tiny_solve(&solver); //연산 
    // unsigned long end = micros();
    // printf("%10.4f %10d %10d\n", (work.x.col(0) - work.Xref.col(0)).squaredNorm(), work.iter, end - start);

    float yaw_moment = work.u.col(0)(0);

    writeToSimulink(yaw_moment);  // 결과를 Simulink로 전송 
    clearSerialBuffer();  //clearbuffer 반드시 있어야함. 
    delay(10);  //delay 반드시 있어야함. 
}
}


