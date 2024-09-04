#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <eigen3/Eigen/Dense>
using namespace std;

#define PI 3.14159265358979323846


struct DHParameter {
    double a;
    double alpha;
    double d;
    double theta;
};

struct Cart_Pose {
    double x;
    double y;
    double z;
    double roll;
    double pitch;
    double yaw;
};

struct Cart_Vel {
    double x_dot;
    double y_dot;
    double z_dot;
    double roll_dot;
    double pitch_dot;
    double yaw_dot;
};

struct Joint_Pose
{
    double Theta_1;
    double Theta_2;
    double Theta_3;
    double Theta_4;
    double Theta_5;
    double Theta_6;

};

struct Joint_Vel
{

double q1_dot;
double q2_dot;
double q3_dot;
double q4_dot;
double q5_dot;
double q6_dot;
};

using Matrix4x4 = std::array<std::array<double, 4>, 4>;
using Vector3 = std::array<double, 3>;
using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

Matrix transposeMatrix(const Matrix& mat) {
    Matrix transposed(mat[0].size(), std::vector<double>(mat.size()));
    for (size_t i = 0; i < mat.size(); ++i) {
        for (size_t j = 0; j < mat[0].size(); ++j) {
            transposed[j][i] = mat[i][j];
        }
    }
    return transposed;
}



Matrix multiplyMatrices(const Matrix& A, const Matrix& B) {
    Matrix result(A.size(), std::vector<double>(B[0].size(), 0.0));
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < B[0].size(); ++j) {
            for (size_t k = 0; k < A[0].size(); ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}


Vector multiplyMatrixVector(const Matrix& mat, const Vector& vec) {
    Vector result(mat.size(), 0.0);
    for (size_t i = 0; i < mat.size(); ++i) {
        for (size_t j = 0; j < vec.size(); ++j) {
            result[i] += mat[i][j] * vec[j];
        }
    }
    return result;
}


Matrix dhTransform(double a, double alpha, double d, double theta) {
    Matrix T =  {{
        {cos(theta), -sin(theta) * cos(alpha),  sin(theta) * sin(alpha), a * cos(theta)},
        {sin(theta),  cos(theta) * cos(alpha), -cos(theta) * sin(alpha), a * sin(theta)},
        {0,           sin(alpha),               cos(alpha),               d},
        {0,           0,                        0,                        1}
    }};
    return T;
}

Vector3 crossProduct(const Vector3& a, const Vector3& b) {
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

Vector3 subtractVectors(const Vector3& a, const Vector3& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Matrix forwardKinematicsTransform(const std::vector<DHParameter>& dhParams) {
    int count=0;
    Matrix T = {{
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    }};
    
    for (const auto& param : dhParams) {
        T = multiplyMatrices(T, dhTransform(param.a, param.alpha, param.d, param.theta));
    }
    return T; 
}

    
Cart_Pose forwardKinematics(const std::vector<DHParameter>& dhParams) {
    Matrix T = forwardKinematicsTransform(dhParams);

    Cart_Pose pose;
    pose.x = T[0][3];
    pose.y = T[1][3];
    pose.z = T[2][3];

    // Extracting roll, pitch, yaw from rotation matrix T
    pose.roll = atan2(T[2][1], T[2][2]);
    pose.pitch = atan2(-T[2][0], sqrt(T[2][1] * T[2][1] + T[2][2] * T[2][2]));//0
    pose.yaw = atan2(T[1][0], T[0][0]);
    return pose;
}

Matrix Inverse_Matrix(Matrix&T){
    Matrix T_inv = {{
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    }};
    Eigen::MatrixXd ans(4,4);
    for(int i =0; i<4;i++){
        for(int j=0;j<4;j++){
            ans(i,j)=T[i][j];
        }
    }
    // Code for pseudo-inverse calculation 
    Eigen::MatrixXd inverse = ans.inverse();
     for(int i =0; i<4;i++){
        for(int j=0;j<4;j++){
            T_inv[i][j]=inverse(i,j);
        }
    }
    return T_inv;
}


std::vector<std::vector<double>> calculateJacobian(const std::vector<DHParameter>& dhParams) {
    std::vector<Matrix> transforms;
    Matrix T = {{
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    }};

    for (const auto& param : dhParams) {
        T = multiplyMatrices(T, dhTransform(param.a, param.alpha, param.d, param.theta));
        transforms.push_back(T);
    }

    Vector3 p_e = {T[0][3], T[1][3], T[2][3]};
    std::vector<std::vector<double>> J(6, std::vector<double>(6, 0.0));

    for (int i = 0; i < 6; ++i) {
        Vector3 z_i = {transforms[i][0][2], transforms[i][1][2], transforms[i][2][2]};
        Vector3 p_i = {transforms[i][0][3], transforms[i][1][3], transforms[i][2][3]};
        Vector3 p_e_minus_p_i = subtractVectors(p_e, p_i);
        Vector3 J_v = crossProduct(z_i, p_e_minus_p_i);

        J[0][i] = J_v[0];
        J[1][i] = J_v[1];
        J[2][i] = J_v[2];

        J[3][i] = z_i[0];
        J[4][i] = z_i[1];
        J[5][i] = z_i[2];
    }

    return J;
}

Matrix pseudoInverse(const Matrix& mat) {
    // pseudo inverse method is computationally time consuming so for real time applications like otg or cartesian jogging it has to be avoided so a new method is required
    Matrix U, V;
    Vector S;
    Matrix transposed = transposeMatrix(mat);
    Matrix ans1(mat.size(),std::vector<double>(mat.size()));
    Matrix ans2(mat.size(),std::vector<double>(mat.size()));
    Matrix result(mat[0].size(), Vector(mat.size()));

    // Pseudo-inverse calculation
   ans1 =  multiplyMatrices(mat,transposed);
    Eigen::MatrixXd ans(ans1.size(),ans1.size());
    for(int i =0; i<ans1.size();i++){
        for(int j=0;j<ans1.size();j++){
            ans(i,j)=ans1[i][j];
        }
    }
    // Code for pseudo-inverse calculation 
    Eigen::MatrixXd inverse = ans.inverse();
     for(int i =0; i<ans2.size();i++){
        for(int j=0;j<ans2.size();j++){
            ans1[i][j]=inverse(i,j);
        }
    }
    result = multiplyMatrices(transposed,ans1);

    return result;
}


void Inverse_position_kinematics(const std::vector<DHParameter>& dhParams, Matrix T){
    //****************/ robot parameters**********
    double d1 = 0.1454;
    double d2 = 0;
    double a2 = 0.147;
    double d3 = 0;
    double a3 = 0.147;
    double d4 = 0.0405+0.0352+0.005;
    double d5 = 0.08;
    double d6 = 0.03;


    double th11 =0;
    double th12 =0;
    double th2 =0;
    double th3 =0;
    double th4 =0;
    double th51 =0;
    double th52 =0;
    double th53 =0;
    double th54 =0;
    double th61 =0;
    double th62 =0;
    double th63 =0;
    double th64 =0;


     Vector3 p6 = {T[0][3],T[1][3],T[2][3]}; // COORDINATES OF EE
    Vector3 z6 = {T[0][2],T[1][2],T[2][2]}; //DIRECTION OF AXIS Z6 OF EE
    Vector3 x6 = {T[0][0],T[1][0],T[2][0]}; //DIRECTION OF AXIS Z6 OF EE
    Vector3 p5 = {T[0][3]-(d6*z6[0]),T[1][3]-(d6*z6[1]),T[2][3]-(d6*z6[2])}; // COORDINATES OF FRAME 5

    // ************Calculating  th1**************
    // atan(val) -> has range (-pi/2,pi/2)
    // ************While atan2(y,x) has range (-pi,pi) and takes into account sign of x,y!******************
    double alpha = atan2(p5[1],p5[0]);
    // double alpha = (atan((p5[1]/p5[0]))); ..
    double beta = (d2+d3+d4)/sqrt(p5[1]*p5[1] + p5[0]*p5[0]);
    std::cout<<"alpha "<<alpha<<" "<<"beta "<<beta<<endl;
    // th11 = alpha - asin(beta);
   th11 = alpha + asin(beta);
   th12 = PI + alpha - asin(beta); // giving correct
//    th12 = atan(p5[1]/p5[0]) - asin(beta);
    
    std::cout<<"th11 "<<th11<<endl;
    std::cout<<"th12 "<<th12<<endl;

    // *********CALCULATING TH5***********
    double j51 = p6[0]*sin(th11) - p6[1]*cos(th11);
    double j52 = p6[0]*sin(th12) - p6[1]*cos(th12);
    double w51 = (j51 - d2 - d3 - d4)/d6;
    double w52 = (j52 - d2 - d3 - d4)/d6;

    // std::cout<<"w51 "<<w51<<endl;
    // std::cout<<"w52 "<<w52<<endl;    
    th51 = acos(w51);
    th52 = -acos(w51);
    th53 = acos(w52);
    th54 = -acos(w52);


    std::cout<<"th51 "<<th51<<endl;
    std::cout<<"th52 "<<th52<<endl;
    std::cout<<"th53 "<<th53<<endl;
    std::cout<<"th54 "<<th54<<endl;

  Matrix  T0_1 = dhTransform(0,PI/2,0.1454,th11);
  Matrix  T0_2 = dhTransform(0,PI/2,0.1454,th12);

  Matrix T1_01 = Inverse_Matrix(T0_1);
  Matrix T1_02 = Inverse_Matrix(T0_2);

  Matrix T1_61 = multiplyMatrices(T1_01,T);
  Matrix T1_62 = multiplyMatrices(T1_02,T);

  Matrix T6_11 = Inverse_Matrix(T1_61);
  Matrix T6_12 = Inverse_Matrix(T1_62);

  double p611 = T6_11[0][2]/sin(th51); //Zx
  double p612 = T6_11[0][2]/sin(th52);
  double p621 = T6_12[0][2]/sin(th53);
  double p622 = T6_12[0][2]/sin(th54);

  double q611 = T6_11[1][2]/sin(th51); //Zy
  double q612 = T6_11[1][2]/sin(th52);
  double q621 = T6_12[1][2]/sin(th53);
  double q622 = T6_12[1][2]/sin(th54);

  th61 = atan2(-q611,p611);
  th62 =  atan2(-q612,p612);
  th63 = atan2(-q621,p621);
  th64 = atan2(-q622,p622);

    std::cout<<"th61 "<<th61<<endl;
    std::cout<<"th62 "<<th62<<endl;
    std::cout<<"th63 "<<th63<<endl;
    std::cout<<"th64 "<<th64<<endl;

    Matrix T4_51 = dhTransform(0, -PI/2, 0.08, th51);
    Matrix T4_52 = dhTransform(0, -PI/2, 0.08, th52);
    Matrix T4_53 = dhTransform(0, -PI/2, 0.08, th53);
    Matrix T4_54 = dhTransform(0, -PI/2, 0.08, th54);

    Matrix T5_61 = dhTransform(0, 0, 0.03, th61);
    Matrix T5_62 = dhTransform(0, 0, 0.03, th62);
    Matrix T5_63 = dhTransform(0, 0, 0.03, th63);
    Matrix T5_64 = dhTransform(0, 0, 0.03, th64);

    Matrix T4_61 = multiplyMatrices(T4_51, T5_61);
    Matrix T4_62 = multiplyMatrices(T4_52, T5_62);
    Matrix T4_63 = multiplyMatrices(T4_53, T5_63);
    Matrix T4_64 = multiplyMatrices(T4_54, T5_64);

    Matrix T6_41 = Inverse_Matrix(T4_61);
    Matrix T6_42 = Inverse_Matrix(T4_62);
    Matrix T6_43 = Inverse_Matrix(T4_63);
    Matrix T6_44 = Inverse_Matrix(T4_64);

    Matrix T1_41 = multiplyMatrices(T1_61, T6_41);
    Matrix T1_42 = multiplyMatrices(T1_61, T6_42);
    Matrix T1_43 = multiplyMatrices(T1_62, T6_43);
    Matrix T1_44 = multiplyMatrices(T1_61, T6_44);

    Vector3 p31 = {T1_41[0][3] - d4*T1_41[0][1], T1_41[1][3] - d4*T1_41[1][1], T1_41[2][3] - d4*T1_41[2][1]};
    Vector3 p32 = {T1_42[0][3] - d4*T1_42[0][1], T1_42[1][3] - d4*T1_42[1][1], T1_42[2][3] - d4*T1_42[2][1]};
    Vector3 p33 = {T1_43[0][3] - d4*T1_43[0][1], T1_43[1][3] - d4*T1_43[1][1], T1_43[2][3] - d4*T1_43[2][1]};
    Vector3 p34 = {T1_44[0][3] - d4*T1_44[0][1], T1_44[1][3] - d4*T1_44[1][1], T1_44[2][3] - d4*T1_44[2][1]};


    double val = 2*a2*a3 + a2*a2 + a3*a3;
    
    if((p31[0]*p31[0] + p31[1]*p31[1] + p31[2]*p31[2])<val || (p32[0]*p32[0] + p32[1]*p32[1] + p32[2]*p32[2])<val || (p33[0]*p33[0] + p33[1]*p33[1] + p33[2]*p33[2])<val || (p34[0]*p34[0] + p34[1]*p34[1] + p34[2]*p34[2])<val)
    {
    double j31 = (-(a2*a2 + a3*a3) + (p31[0]*p31[0] + p31[1]*p31[1] + p31[2]*p31[2])) / (2*(a2*a3));
    double j32 = (-(a2*a2 + a3*a3) + (p32[0]*p32[0] + p32[1]*p32[1] + p32[2]*p32[2])) / (2*(a2*a3));
    double j33 = (-(a2*a2 + a3*a3) + (p33[0]*p33[0] + p33[1]*p33[1] + p33[2]*p33[2])) / (2*(a2*a3));
    double j34 = (-(a2*a2 + a3*a3) + (p34[0]*p34[0] + p34[1]*p34[1] + p34[2]*p34[2])) / (2*(a2*a3));

    double th31 = acos(j31);
    double th32 = acos(j32);
    double th33 = acos(j33);
    double th34 = acos(j34);
    double th35 = -acos(j31);
    double th36 = -acos(j32);
    double th37 = -acos(j33);
    double th38 = -acos(j34);
    cout<<endl;
    cout << "th31 " << th31 << endl;
    cout << "th32 " << th32 << endl;
    cout << "th33 " << th33 << endl;
    cout << "th34 " << th34 << endl;
    cout << "th35 " << th35 << endl;
    cout << "th36 " << th36 << endl;
    cout << "th37 " << th37 << endl;
    cout << "th38 " << th38 << endl;

    double th21 = atan2(p31[1], p31[0]) - asin((a3 * sin(th31)) / sqrt(p31[0]*p31[0] + p31[1]*p31[1] + p31[2]*p31[2]));
    double th22 = atan2(p32[1], p32[0]) - asin((a3 * sin(th32)) / sqrt(p32[0]*p32[0] + p32[1]*p32[1] + p32[2]*p32[2]));
    double th23 = atan2(p33[1], p33[0]) - asin((a3 * sin(th33)) / sqrt(p33[0]*p33[0] + p33[1]*p33[1] + p33[2]*p33[2]));
    double th24 = atan2(p34[1], p34[0]) - asin((a3 * sin(th34)) / sqrt(p34[0]*p34[0] + p34[1]*p34[1] + p34[2]*p34[2]));
    double th25 = atan2(p31[1], p31[0]) - asin((a3 * sin(th35)) / sqrt(p31[0]*p31[0] + p31[1]*p31[1] + p31[2]*p31[2]));
    double th26 = atan2(p32[1], p32[0]) - asin((a3 * sin(th36)) / sqrt(p32[0]*p32[0] + p32[1]*p32[1] + p32[2]*p32[2]));
    double th27 = atan2(p33[1], p33[0]) - asin((a3 * sin(th33)) / sqrt(p33[0]*p33[0] + p33[1]*p33[1] + p33[2]*p33[2]));
    double th28 = atan2(p34[1], p34[0]) - asin((a3 * sin(th34)) / sqrt(p34[0]*p34[0] + p34[1]*p34[1] + p34[2]*p34[2]));

    if(th31<-2*PI){
        th21 = 2*PI + th21;
    }

    if(th32<-2*PI){
        th22 = 2*PI + th22;
    }

    if(th33<-2*PI){
        th23 = 2*PI + th23;
    }

    if(th34<-2*PI){
        th24 = 2*PI + th24;
    }
    if(th35<-2*PI){
        th25 = 2*PI + th25;
    }
    if(th36<-2*PI){
        th26 = 2*PI + th26;
    }
    if(th37<-2*PI){
        th27 = 2*PI + th27;
    }
    if(th38<-2*PI){
        th28 = 2*PI + th28;
    }

    cout<<endl;
    cout << "th21 " << th21 << endl;
    cout << "th22 " << th22 << endl;
    cout << "th23 " << th23 << endl;
    cout << "th24 " << th24 << endl;
    cout << "th25 " << th25 << endl;
    cout << "th26 " << th26 << endl;
    cout << "th27 " << th27 << endl;
    cout << "th28 " << th28 << endl;

    Matrix T1_21 = dhTransform(0.147, 0, 0, th21);
    Matrix T1_22 = dhTransform(0.147, 0, 0, th22);
    Matrix T1_23 = dhTransform(0.147, 0, 0, th23);
    Matrix T1_24 = dhTransform(0.147, 0, 0, th24);
    Matrix T1_25 = dhTransform(0.147, 0, 0, th25);
    Matrix T1_26 = dhTransform(0.147, 0, 0, th26);
    Matrix T1_27 = dhTransform(0.147, 0, 0, th27);
    Matrix T1_28 = dhTransform(0.147, 0, 0, th28);

    Matrix T2_31 = dhTransform(0.147, 0, 0, th31);
    Matrix T2_32 = dhTransform(0.147, 0, 0, th32);
    Matrix T2_33 = dhTransform(0.147, 0, 0, th33);
    Matrix T2_34 = dhTransform(0.147, 0, 0, th34);
    Matrix T2_35 = dhTransform(0.147, 0, 0, th35);
    Matrix T2_36 = dhTransform(0.147, 0, 0, th36);
    Matrix T2_37 = dhTransform(0.147, 0, 0, th37);
    Matrix T2_38 = dhTransform(0.147, 0, 0, th38);

    Matrix T1_31 = multiplyMatrices(T1_21, T2_31);
    Matrix T1_32 = multiplyMatrices(T1_22, T2_32);
    Matrix T1_33 = multiplyMatrices(T1_23, T2_33);
    Matrix T1_34 = multiplyMatrices(T1_24, T2_34);
    Matrix T1_35 = multiplyMatrices(T1_25, T2_35);
    Matrix T1_36 = multiplyMatrices(T1_26, T2_36);
    Matrix T1_37 = multiplyMatrices(T1_27, T2_37);
    Matrix T1_38 = multiplyMatrices(T1_28, T2_38);

    Matrix T3_11 = Inverse_Matrix(T1_31);
    Matrix T3_12 = Inverse_Matrix(T1_32);
    Matrix T3_13 = Inverse_Matrix(T1_33);
    Matrix T3_14 = Inverse_Matrix(T1_34);
    Matrix T3_15 = Inverse_Matrix(T1_35);
    Matrix T3_16 = Inverse_Matrix(T1_36);
    Matrix T3_17 = Inverse_Matrix(T1_37);
    Matrix T3_18 = Inverse_Matrix(T1_38);

    Matrix T3_41 = multiplyMatrices(T3_11,T1_41);
    Matrix T3_42 = multiplyMatrices(T3_12,T1_42);
    Matrix T3_43 = multiplyMatrices(T3_13,T1_43);
    Matrix T3_44 = multiplyMatrices(T3_14,T1_44);
    Matrix T3_45 = multiplyMatrices(T3_15,T1_41);
    Matrix T3_46 = multiplyMatrices(T3_16,T1_42);
    Matrix T3_47 = multiplyMatrices(T3_17,T1_43);
    Matrix T3_48 = multiplyMatrices(T3_18,T1_44);

    double th41 = atan2(T3_41[1][0], T3_41[0][0]);
    double th42 = atan2(T3_42[1][0], T3_42[0][0]);
    double th43 = atan2(T3_43[1][0], T3_43[0][0]);
    double th44 = atan2(T3_44[1][0], T3_44[0][0]);
    double th45 = atan2(T3_45[1][0], T3_45[0][0]);
    double th46 = atan2(T3_46[1][0], T3_46[0][0]);
    double th47 = atan2(T3_47[1][0], T3_47[0][0]);
    double th48 = atan2(T3_48[1][0], T3_48[0][0]);

    cout<<endl;
    cout << "th41 " << th41 << endl;
    cout << "th42 " << th42 << endl;
    cout << "th43 " << th43 << endl;
    cout << "th44 " << th44 << endl;
    cout << "th45 " << th45 << endl;
    cout << "th46 " << th46 << endl;
    cout << "th47 " << th47 << endl;
    cout << "th48 " << th48 << endl;


    }
    else{
        cout<<"***************no sol***************"<<endl;
    }

    

}


int main() {
    Joint_Pose jpose;
    jpose = {0.56,1.24,1.34,-0.9,0.5,1.59};
    std::vector<DHParameter> dhParams = {
        {0,     PI / 2,    0.1454,    jpose.Theta_1},
        {0.147,    0,         0,      jpose.Theta_2},
        {0.147,   0,    0,   jpose.Theta_3},
        {0,      PI / 2,    0.0405+0.0352+0.005,    jpose.Theta_4},
        {0,    -PI / 2,    0.08,      jpose.Theta_5},
        {0,    0,         0.03,      jpose.Theta_6}
    };

    Cart_Pose pose = forwardKinematics(dhParams);
    Matrix T = forwardKinematicsTransform(dhParams);
    std::cout << "Position: (" << pose.x << ", " << pose.y << ", " << pose.z << ")\n";
    // std::cout << "Orientation: (roll: " << pose.roll << ", pitch: " << pose.pitch << ", yaw: " << pose.yaw << ")\n";
    Inverse_position_kinematics(dhParams,T);
    
    return 0;
}
