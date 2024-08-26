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




void Inverse_position_kinematics(const std::vector<DHParameter>& dhParams, Matrix T){
    //****************/ robot parameters**********
    double d1 = 0.1454;
    double d2 = 0.0352;
    double a2 = 0.147;
    double d3 = 0.005;
    double a3 = 0.147;
    double d4 = 0.0405;
    double d5 = 0.08;
    double d6 = 0.03;

    // ************joint parameters***************
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
    double th65 =0;
    double th66 =0;
    double th67 =0;
    double th68 =0;

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

    //***************CALCULATING TH6********************
    double a61 = cos(th51)*sin(th11);
    double a62 = cos(th52)*sin(th11); // no meaning
    double a63 = cos(th53)*sin(th12);
    double a64 = cos(th54)*sin(th12); // no meaning

    double b61 = cos(th11)*sin(th51);
    double b62 = cos(th11)*sin(th52);
    double b63 = cos(th12)*sin(th53);
    double b64 = cos(th12)*sin(th54);

    double p61 = (a61 - z6[0])/b61;
    double p62 = (a62 - z6[0])/b62;
    double p63 = (a63 - z6[0])/b63;
    double p64 = (a64 - z6[0])/b64;

    double q61 = -z6[2]/sin(th51);
    double q62 = -z6[2]/sin(th52);
    double q63 = -z6[2]/sin(th53);
    double q64 = -z6[2]/sin(th54);

    double r61 = q61*cos(th51);
    double r62 = q62*cos(th52);
    double r63 = q63*cos(th53);
    double r64 = q64*cos(th54);

    double y61 = acos(p61/sqrt(p61*p61 + r61*r61));
    double y62 = acos(p62/sqrt(p62*p62 + r62*r62));
    double y63 = acos(p63/sqrt(p63*p63 + r63*r63));
    double y64 = acos(p64/sqrt(p64*p64 + r64*r64));

    double m61 = x6[2]/sqrt(p61*p61 + r61*r61);
    double m62 = x6[2]/sqrt(p62*p62 + r62*r62);
    double m63 = x6[2]/sqrt(p63*p63 + r63*r63);
    double m64 = x6[2]/sqrt(p64*p64 + r64*r64);

    // th62 = -PI + asin(m62) - y62;
    // th61 = -PI + asin(m61) - y61;
    // th63 = -PI + asin(m63) - y63;
    // th64 = -PI + asin(m64) - y64;

    th61 = asin(m61) - y61;
    th62 = asin(m62) - y62;
    th63 = asin(m63) - y63;
    th64 = asin(m64) - y64;
    th65 = asin(m61) + y61;
    th66 = asin(m62) + y62;
    th67 = asin(m63) + y63;
    th68 = asin(m64) + y64;
    // th69 = asin(m61) - y61;
    // th610 = asin(m62) - y62;
    // th611 = asin(m63) - y63;
    // th612 = asin(m64) - y64;
    // th613 = asin(m61) + y61;
    // th614 = asin(m62) + y62;
    // th615 = asin(m63) + y63;
    // th616= asin(m64) + y64;

    std::cout<<"th61 "<<th61<<endl;
    std::cout<<"th62 "<<th62<<endl;
    std::cout<<"th63 "<<th63<<endl;
    std::cout<<"th64 "<<th64<<endl;
    std::cout<<"th65 "<<th65<<endl;
    std::cout<<"th66 "<<th66<<endl;
    std::cout<<"th67 "<<th67<<endl;
    std::cout<<"th68 "<<th68<<endl;



}


int main() {
    Joint_Pose jpose;
    jpose = {1.02,3,1.56,0.67,-1.56,1.44};
    std::vector<DHParameter> dhParams = {
        {0,     PI / 2,    0.1454,    jpose.Theta_1},
        {0.147,    0,         0.0352,      jpose.Theta_2},
        {0.147,   0,    0.005,   jpose.Theta_3},
        {0,      PI / 2,    0.0405,    jpose.Theta_4},
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