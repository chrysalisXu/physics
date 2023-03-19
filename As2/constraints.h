#ifndef CONSTRAINTS_HEADER_FILE
#define CONSTRAINTS_HEADER_FILE

using namespace Eigen;
using namespace std;

typedef enum ConstraintType{DISTANCE, COLLISION} ConstraintType;   //You can expand it for more constraints
typedef enum ConstraintEqualityType{EQUALITY, INEQUALITY} ConstraintEqualityType;

//there is such constraints per two variables that are equal. That is, for every attached vertex there are three such constraints for (x,y,z);
class Constraint{
public:
  
  int m1, m2;                     //Two participating meshes (can be the same)  - auxiliary data for users (constraint class shouldn't use that)
  int v1, v2;                     //Two vertices from the respective meshes - auxiliary data for users (constraint class shouldn't use that)
  double invMass1, invMass2;       //inverse masses of two bodies
  double refValue;                //Reference values to use in the constraint, when needed (like distance), d12
  RowVector3d refVector;             //Reference vector when needed (like vector) , n hat
  double CRCoeff;                 //extra velocity bias
  ConstraintType constraintType;  //The type of the constraint, and will affect the value and the gradient. This SHOULD NOT change after initialization!
  ConstraintEqualityType constraintEqualityType;  //whether the constraint is an equality or an inequality
  
  Constraint(const ConstraintType _constraintType, const ConstraintEqualityType _constraintEqualityType, const int& _m1, const int& _v1, const int& _m2, const int& _v2, const double& _invMass1, const double& _invMass2, const RowVector3d& _refVector, const double& _refValue, const double& _CRCoeff):constraintType(_constraintType), constraintEqualityType(_constraintEqualityType), m1(_m1), v1(_v1), m2(_m2), v2(_v2), invMass1(_invMass1), invMass2(_invMass2),  refValue(_refValue), CRCoeff(_CRCoeff){
    refVector=_refVector;
  }
  
  ~Constraint(){}
  
  
  
  //computes the impulse needed for all particles to resolve the velocity constraint, and corrects the velocities accordingly.
  //The velocities are a vector (vCOM1, w1, vCOM2, w2) in both input and output.
  //returns true if constraint was already valid with "currVelocities", and false otherwise (false means there was a correction done)
  //currCOMPositions is a 2x3 matrix, where each row is per one of the sides of the constraints; the rest of the relevant variables are similar, and so should the outputs be resized.
  bool resolveVelocityConstraint(const MatrixXd& currCOMPositions, const MatrixXd& currVertexPositions, const MatrixXd& currCOMVelocities, const MatrixXd& currAngularVelocities, const Matrix3d& invInertiaTensor1, const Matrix3d& invInertiaTensor2, MatrixXd& correctedCOMVelocities, MatrixXd& correctedAngularVelocities, double tolerance){
    
    MatrixXd invMassMatrix=MatrixXd::Zero(12,12);
    RowVectorXd constGradient(12); // constraint Gradient is Jacobian 1*12
    
    
    /**************
     TODO: write velocity correction procedure:
     1. If the velocity Constraint is satisfied up to tolerate ("abs(Jv)<=tolerance"), set corrected values to original ones and return true
     
     2. Otherwise, correct linear and angular velocities as learnt in class.
     
     Note to differentiate between different constraint types; for inequality constraints you don't do anything unless it's unsatisfied.
     ***************/
    Matrix3d IdentityMat3d = Matrix3d::Identity(); // 3*3 identity matrix
    //invMassMatrix = (invMass1 * IdentityMat3d, invInertiaTensor1, invMass2 * IdentityMat3d, invInertiaTensor2);
    invMassMatrix.block(0, 0, 3, 3) = invMass1 * IdentityMat3d;
    invMassMatrix.block(3, 3, 3, 3) = invInertiaTensor1;
    invMassMatrix.block(6, 6, 3, 3) = invMass2 * IdentityMat3d;
    invMassMatrix.block(9, 9, 3, 3) = invInertiaTensor2;
    double Constraint_x1_x2 = (currCOMPositions.row(0) - currCOMPositions.row(1)).norm() - refValue;
    RowVectorXd velocity_Vector(12); // 1*12
    velocity_Vector = (currCOMVelocities.row(0), currAngularVelocities.row(0), currCOMVelocities.row(1), currAngularVelocities.row(1));
    RowVector3d n_hat = (currCOMPositions.row(0) - currCOMPositions.row(1))/(currCOMPositions.row(0) - currCOMPositions.row(1)).norm();
    RowVector3d rA = currCOMPositions.row(0) - currVertexPositions.row(0);
    RowVector3d rB = currCOMPositions.row(1) - currVertexPositions.row(1);
    constGradient = (n_hat, rA.cross(n_hat), -n_hat, rB.cross(n_hat)); // 1*12
    // define lamda: Lagrange multiplier, lamda_up: numerator, lamda_down: denominator
    double lamda_up = constGradient.transpose().dot(velocity_Vector); // 1*12 .dot 1*12
    double lamda_down = constGradient * invMassMatrix * constGradient.transpose(); // 1*12 * 12*12 * 12*1
    double lamda = -lamda_up / lamda_down;
    RowVectorXd delta_Velocity(12);
    delta_Velocity = -lamda * invMassMatrix * constGradient.transpose(); //12*12 * 12*1

    // RowVector3d impulse = constGradient.transpose() * lamda;

    // need to determine Δv such that Ji · (v + Δv) = 0
    //if (constGradient.dot(velocity_Vector + delta_Velocity) == 0) {
    if (abs(constGradient.dot(velocity_Vector)) <= tolerance) {
        correctedCOMVelocities=currCOMVelocities;
        correctedAngularVelocities=currAngularVelocities;
        return true;
    }
    else {
        RowVectorXd correctedVelocity_Vector = velocity_Vector + delta_Velocity;
        correctedCOMVelocities = correctedVelocity_Vector.segment(0,3) + correctedVelocity_Vector.segment(6,9);
        correctedAngularVelocities = correctedVelocity_Vector.segment(3, 6) + correctedVelocity_Vector.segment(9, 12);
        return false;
    };
    
  }
  
  //projects the position unto the constraint
  //returns true if constraint was already valid with "currPositions"
  bool resolvePositionConstraint(const MatrixXd& currCOMPositions, const MatrixXd& currConstPositions, MatrixXd& correctedCOMPositions, double tolerance){
    
    MatrixXd invMassMatrix=MatrixXd::Zero(6,6);
    
    
    /**************
     TODO: write position correction procedure:
     1. If the position Constraint is satisfied up to tolerate ("abs(C(p)<=tolerance"), set corrected values to original ones and return true
     
     2. Otherwise, correct COM position as learnt in class. Note that since this is a linear correction, correcting COM position == correcting all positions the same offset. the currConstPositions are used to measure the constraint, and the COM values are corrected accordingly to create the effect.
     
     Note to differentiate between different constraint types; for inequality constraints you don't do anything unless it's unsatisfied.
     ***************/

    //In Position correction, ignore angular velocity
    // So J is also a 1*6 vector
    Matrix3d IdentityMat3d = Matrix3d::Identity(); // 3*3 identity matrix
    invMassMatrix.block(0, 0, 3, 3) = invMass1 * IdentityMat3d;
    invMassMatrix.block(3, 3, 3, 3) = invMass2 * IdentityMat3d;
    
    double Constraint_x1_x2 = (currCOMPositions.row(0) - currCOMPositions.row(1)).norm() - refValue;
    // equality constraint
    if (abs(Constraint_x1_x2) <= tolerance) { 
        correctedCOMPositions = currCOMPositions;
        return true;
    }
    else {
        // correct COM positions
        // according to topic 6, slide 9
        RowVector3d n_hat = (currCOMPositions.row(0) - currCOMPositions.row(1)) / (currCOMPositions.row(0) - currCOMPositions.row(1)).norm();
        RowVector3d r1 = currCOMPositions.row(0);
        RowVector3d r2 = currCOMPositions.row(1);
        RowVector3d delta_r1 = - invMass2/(invMass1 + invMass2) * ((r1 - r2).norm() - refValue)* n_hat;
        RowVector3d delta_r2 = invMass1 / (invMass1 + invMass2) * ((r1 - r2).norm() - refValue) * n_hat;
        correctedCOMPositions.row(0) = r1 + delta_r1;
        correctedCOMPositions.row(1) = r2 + delta_r2;

        return false;
    }
    
    

  }
};



#endif /* constraints_h */
