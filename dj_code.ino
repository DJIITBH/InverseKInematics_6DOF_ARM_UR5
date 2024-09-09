#include <Arduino.h>
#include <math.h>
#include <vector>

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

struct Joint_Pose {
    double Theta_1;
    double Theta_2;
    double Theta_3;
    double Theta_4;
    double Theta_5;
    double Theta_6;
};

struct Joint_Vel {
    double q1_dot;
    double q2_dot;
    double q3_dot;
    double q4_dot;
    double q5_dot;
    double q6_dot;
};

using Matrix = std::vector<std::vector<double>>;
using Vector3 = std::array<double, 3>;
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
    Matrix T = {
        {cos(theta), -sin(theta) * cos(alpha),  sin(theta) * sin(alpha), a * cos(theta)},
        {sin(theta),  cos(theta) * cos(alpha), -cos(theta) * sin(alpha), a * sin(theta)},
        {0,           sin(alpha),               cos(alpha),               d},
        {0,           0,                        0,                        1}
    };
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
    Matrix T = {
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
    
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

Matrix getMinor(const Matrix& mat, int row, int col) {
    Matrix minor(3, std::vector<double>(3));
    int mi = 0, mj = 0;

    for (int i = 0; i < 4; ++i) {
        if (i == row) continue;
        mj = 0;
        for (int j = 0; j < 4; ++j) {
            if (j == col) continue;
            minor[mi][mj] = mat[i][j];
            ++mj;
        }
        ++mi;
    }

    return minor;
}

// Helper function to calculate the determinant of a 3x3 matrix
double determinant3x3(const Matrix& mat) {
    return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
         - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
         + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
}

// Function to calculate the determinant of a 4x4 matrix
double determinant4x4(const Matrix& mat) {
    double det = 0.0;
    for (int i = 0; i < 4; ++i) {
        Matrix minor = getMinor(mat, 0, i);
        double sign = (i % 2 == 0) ? 1.0 : -1.0;
        det += sign * mat[0][i] * determinant3x3(minor);
    }
    return det;
}

// Function to calculate the adjugate of a 4x4 matrix
Matrix adjugate(const Matrix& mat) {
    Matrix adj(4, std::vector<double>(4));
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            Matrix minor = getMinor(mat, i, j);
            double sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;
            adj[j][i] = sign * determinant3x3(minor);
        }
    }
    return adj;
}

// Function to compute the inverse of a 4x4 matrix
Matrix inverse(const Matrix& mat) {
    double det = determinant4x4(mat);
    if (det == 0) {
        throw std::runtime_error("Matrix is not invertible");
    }

    Matrix adj = adjugate(mat);
    Matrix inv(4, std::vector<double>(4));

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            inv[i][j] = adj[i][j] / det;
        }
    }

    return inv;
}

Cart_Pose fk(const Matrix& T) {
    Cart_Pose pose;
    pose.x = T[0][3];
    pose.y = T[1][3];
    pose.z = T[2][3];

    // Extracting roll, pitch, yaw from rotation matrix T
    pose.roll = atan2(T[2][1], T[2][2]);
    pose.pitch = atan2(-T[2][0], sqrt(T[2][1] * T[2][1] + T[2][2] * T[2][2]));
    pose.yaw = atan2(T[1][0], T[0][0]);
    return pose;
}

double esp(const Joint_Pose& jpose, const Cart_Pose& target) {
    std::vector<DHParameter> dhParams = {
        {0,     PI / 2,    0.1454,    jpose.Theta_1},
        {0.147,    0,         0,      jpose.Theta_2},
        {0.147,   0,    0,   jpose.Theta_3},
        {0,      PI / 2,    0.0405 + 0.0352 + 0.005,    jpose.Theta_4},
        {0,    -PI / 2,    0.08,      jpose.Theta_5},
        {0,    0,         0.03,      jpose.Theta_6}
    };
    Cart_Pose pose = forwardKinematics(dhParams);
    double esp[6] = {0, 0, 0, 0, 0, 0};
    esp[0] = abs(pose.x - target.x);
    esp[1] = abs(pose.y - target.y);
    esp[2] = abs(pose.z - target.z);
    esp[3] = abs(pose.roll - target.roll);
    esp[4] = abs(pose.pitch - target.pitch);
    esp[5] = abs(pose.yaw - target.yaw);

    double max = esp[0];
    for (int i = 1; i < 6; i++) {
        if (esp[i] > max) {
            max = esp[i];
        }
    }
    return max; //returns max error! that is between given pose and target pose!
}

Joint_Pose Inverse_Position_Kinematics(const  Matrix &T) {
    // Robot parameters
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
    // std::cout<<"alpha "<<alpha<<" "<<"beta "<<beta<<endl;
    // Serial.print
    // th11 = alpha - asin(beta);
   th11 = alpha + asin(beta);
   th12 = PI + alpha - asin(beta); // giving correct
//    th12 = atan(p5[1]/p5[0]) - asin(beta);
    

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
    // Serial.println(th51);
    // Serial.println(th52);
    // Serial.println(th53);
    // Serial.println(th54);


  Matrix  T0_1 = dhTransform(0,PI/2,0.1454,th11);
  Matrix  T0_2 = dhTransform(0,PI/2,0.1454,th12);

  Matrix T1_01 = inverse(T0_1);
  Matrix T1_02 = inverse(T0_2);

  Matrix T1_61 = multiplyMatrices(T1_01,T);
  Matrix T1_62 = multiplyMatrices(T1_02,T);

  Matrix T6_11 = inverse(T1_61);
  Matrix T6_12 = inverse(T1_62);

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

  // Serial.println(th61);
  //   Serial.println(th62);
  //   Serial.println(th63);
  //   Serial.println(th64);



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

    Matrix T6_41 = inverse(T4_61);
    Matrix T6_42 = inverse(T4_62);
    Matrix T6_43 = inverse(T4_63);
    Matrix T6_44 = inverse(T4_64);

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

    // Serial.println(th31);
    // Serial.println(th32);
    // Serial.println(th33);
    // Serial.println(th34);
    // Serial.println(th35);
    // Serial.println(th36);
    // Serial.println(th37);
    // Serial.println(th38);

    

    double th21 = atan2(p31[1], p31[0]) - asin((a3 * sin(th31)) / sqrt(p31[0]*p31[0] + p31[1]*p31[1] + p31[2]*p31[2]));
    double th22 = atan2(p32[1], p32[0]) - asin((a3 * sin(th32)) / sqrt(p32[0]*p32[0] + p32[1]*p32[1] + p32[2]*p32[2]));
    double th23 = atan2(p33[1], p33[0]) - asin((a3 * sin(th33)) / sqrt(p33[0]*p33[0] + p33[1]*p33[1] + p33[2]*p33[2]));
    double th24 = atan2(p34[1], p34[0]) - asin((a3 * sin(th34)) / sqrt(p34[0]*p34[0] + p34[1]*p34[1] + p34[2]*p34[2]));
    double th25 = atan2(p31[1], p31[0]) - asin((a3 * sin(th35)) / sqrt(p31[0]*p31[0] + p31[1]*p31[1] + p31[2]*p31[2]));
    double th26 = atan2(p32[1], p32[0]) - asin((a3 * sin(th36)) / sqrt(p32[0]*p32[0] + p32[1]*p32[1] + p32[2]*p32[2]));
    double th27 = atan2(p33[1], p33[0]) - asin((a3 * sin(th33)) / sqrt(p33[0]*p33[0] + p33[1]*p33[1] + p33[2]*p33[2]));
    double th28 = atan2(p34[1], p34[0]) - asin((a3 * sin(th34)) / sqrt(p34[0]*p34[0] + p34[1]*p34[1] + p34[2]*p34[2]));

    // Serial.println(th21);
    // Serial.println(th22);
    // Serial.println(th23);
    // Serial.println(th24);
    // Serial.println(th25);
    // Serial.println(th26);
    // Serial.println(th27);
    // Serial.println(th28);


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

    Matrix T3_11 = inverse(T1_31);
    Matrix T3_12 = inverse(T1_32);
    Matrix T3_13 = inverse(T1_33);
    Matrix T3_14 = inverse(T1_34);
    Matrix T3_15 = inverse(T1_35);
    Matrix T3_16 = inverse(T1_36);
    Matrix T3_17 = inverse(T1_37);
    Matrix T3_18 = inverse(T1_38);

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


// ****************THE 8 SOLUTION SET*******************
    // Joint_Pose jpose1 = {th11,th51,th61,th31,th21,th41};
    Joint_Pose jpose1 = {th11,th21,th31,th41,th51,th61};
    // Joint_Pose jpose2 = {th11,th51,th61,th35,th25,th45};
    Joint_Pose jpose2 = {th11,th25,th35,th45,th51,th61};
    // Joint_Pose jpose3 = {th11,th52,th62,th32,th22,th42};
    Joint_Pose jpose3 = {th11,th22,th32,th42,th52,th62};
    // Joint_Pose jpose4 = {th11,th52,th62,th36,th26,th46};
    Joint_Pose jpose4 = {th11,th26,th36,th46,th52,th62};
    // Joint_Pose jpose5 = {th12,th53,th63,th33,th23,th43};
    Joint_Pose jpose5 = {th12,th23,th33,th43,th53,th63};
    // Joint_Pose jpose6 = {th12,th53,th63,th37,th27,th47};
    Joint_Pose jpose6 = {th12,th27,th37,th47,th53,th63};
    // Joint_Pose jpose7 = {th12,th54,th64,th34,th24,th44};
    Joint_Pose jpose7 = {th12,th24,th34,th44,th54,th64};
    // Joint_Pose jpose8 = {th12,th54,th64,th38,th28,th48};
    Joint_Pose jpose8 = {th12,th28,th38,th48,th54,th64};

    Cart_Pose pose_target = fk(T);

    // Serial.println(jpose2.Theta_2);
    // Serial.println(jpose3.Theta_6);
    // Serial.print("Theta1 : "); 
    // Serial.print(jpose1.Theta_1);
    // Serial.print("  Theta2 : ");
    // Serial.print(jpose1.Theta_2);
    // Serial.print("  Theta3 : ");
    // Serial.print(jpose1.Theta_3);
    // Serial.print( "   Theta4 : ");
    // Serial.print(jpose1.Theta_4);
    // Serial.print("  Theta5 : ");
    // Serial.print(jpose1.Theta_5);
    // Serial.print("  Theta6 : ");
    // Serial.println(jpose1.Theta_6);

    // Serial.println(jpose4.Theta_6);
    // Serial.println(jpose5.Theta_6);
    // Serial.println(jpose6.Theta_6);
    // Serial.println(jpose7.Theta_6);
    // Serial.println(jpose8.Theta_6);

    // *************FIND ERRORS OF SOLUTIONS**********
    double esp1 = esp(jpose1,pose_target);
    double esp2 = esp(jpose2,pose_target);
    double esp3 = esp(jpose3,pose_target);
    double esp4 = esp(jpose4,pose_target);
    double esp5 = esp(jpose5,pose_target);
    double esp6 = esp(jpose6,pose_target);
    double esp7 = esp(jpose7,pose_target);
    double esp8 = esp(jpose8,pose_target);

    std::vector<double> esp = {esp1, esp2, esp3, esp4, esp5, esp6, esp7, esp8};
    double min = esp1; //FINDING MINIMUM ERROR FROM THOSE MAX ERRORS! GETTING BEST POSSIBLE SOLUTIONS
    std::pair<double,int> p;
    p = {min,1};
    for(int i=0;i<8;i++){
        if(esp[i]<min){
            min = esp[i];
            p = {min,i+1};
        }
    }
    int c = p.second;
   switch (c) {
        case 1:
            return jpose1;
            break;
        case 2:
             return jpose2;
            break;
        case 3:
             return jpose3;
            break;
        case 4:
             return jpose4;
            break;
        case 5:
            return jpose5;
            break;
        case 6:
            return jpose6;
            break;
        case 7:
             return jpose7;
            break;
        case 8:
             return jpose8;
            break;
        default:
             return jpose1;
            break;
    }

    }
    else{
      Joint_Pose j = {0,PI/2,0,PI/2,0,0};
      return j;
    }
    
    
 // Adjust according to the minimum error
}


void setup(){
  Serial.begin(115200);

}

void loop(){

    Joint_Pose jpose;
    jpose = {0.01, 0.02, 0.03, 0.03, 0.02, 0.01};

    // Define the DH parameters
    std::vector<DHParameter> dhParams = {
        {0, PI / 2, 0.1454, jpose.Theta_1},
        {0.147, 0, 0, jpose.Theta_2},
        {0.147, 0, 0, jpose.Theta_3},
        {0, PI / 2, 0.0405 + 0.0352 + 0.005, jpose.Theta_4},
        {0, -PI / 2, 0.08, jpose.Theta_5},
        {0, 0, 0.03, jpose.Theta_6}
    };

    // Calculate the forward kinematics and transformation matrix
    Cart_Pose pose = forwardKinematics(dhParams);
    Matrix T = forwardKinematicsTransform(dhParams);

    // Print position
    // Serial.print("Position: (");
    // Serial.print(pose.x);
    // Serial.print(", ");
    // Serial.print(pose.y);
    // Serial.print(", ");
    // Serial.print(pose.z);
    // Serial.println(")");

    // Calculate the inverse kinematics
    Joint_Pose j = Inverse_Position_Kinematics(T);

    // Print joint angles
    Serial.print("Theta_1: ");
    Serial.print(j.Theta_1);
    Serial.print(" Theta_2: ");
    Serial.print(j.Theta_2);
    Serial.print(" Theta_3: ");
    Serial.print(j.Theta_3);
    Serial.print(" Theta_4: ");
    Serial.print(j.Theta_4);
    Serial.print(" Theta_5: ");
    Serial.print(j.Theta_5);
    Serial.print(" Theta_6: ");
    Serial.println(j.Theta_6);

    // Add a delay for readability
    delay(1000);  // Adjust the delay as needed
}
// Add similar conversion for the Inverse_position_kinematics function if needed
