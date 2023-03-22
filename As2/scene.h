#ifndef SCENE_HEADER_FILE
#define SCENE_HEADER_FILE

#include <vector>
#include <fstream>
#include <sstream>
#include <igl/bounding_box.h>
#include <igl/readMESH.h>
#include "ccd.h"
#include "volInt.h"
#include "auxfunctions.h"
#include "constraints.h"

using namespace Eigen;
using namespace std;


void support(const void *_obj, const ccd_vec3_t *_d, ccd_vec3_t *_p);
void stub_dir(const void *obj1, const void *obj2, ccd_vec3_t *dir);
void center(const void *_obj, ccd_vec3_t *dir);



//Impulse is defined as a pair <position, direction>
typedef std::pair<RowVector3d,RowVector3d> Impulse;


// < 10 ^ -5 will be defined as zero
double ZERO_PRECISION = 0.00001;
void precisionCheckZero(Vector3d& target){
  for (int i=0; i<3; i++)
    if (abs(target(i,0)) <= ZERO_PRECISION) 
      target(i,0) = 0;
}


//the class the contains each individual rigid objects and their functionality
class Mesh{
public:
  MatrixXd origV;   //original vertex positions, where COM=(0.0,0.0,0.0) - never change this!
  MatrixXd currV;   //current vertex position
  MatrixXi F;   //faces of the tet mesh
  MatrixXi T;   //Tets in the tet mesh
  
  VectorXi boundTets;  //indices (from T) of just the boundary tets, for collision
  
  //position of object in space. We must always have that currV = QRot(origV, orientation)+ COM
  RowVector4d orientation; //current orientation
  RowVector3d COM;  //current center of mass
  Matrix3d invIT;  //Original *inverse* inertia tensor around the COM, defined in the rest state to the object (so to the canonical world system)
  
  VectorXd tetVolumes;    //|T|x1 tetrahedra volumes
  VectorXd invMasses;     //|T|x1 tetrahedra *inverse* masses
  
  //kinematics
  bool isFixed;  //is the object immobile
  double totalMass;  //sum(1/invMass)
  double totalVolume;
  RowVector3d comVelocity;  //the linear velocity of the center of mass
  RowVector3d angVelocity;  //the angular velocity of the object.
  
  //dynamics
  std::vector<Impulse> currImpulses;  //current list of impulses, updated by collision handling

  //checking collision between bounding boxes, and consequently the boundary tets if succeeds.
  //you do not need to update these functions (isBoxCollide and isCollide) unless you are doing a different collision
  
  bool isBoxCollide(const Mesh& m){
    RowVector3d VMin1=currV.colwise().minCoeff();
    RowVector3d VMax1=currV.colwise().maxCoeff();
    RowVector3d VMin2=m.currV.colwise().minCoeff();
    RowVector3d VMax2=m.currV.colwise().maxCoeff();
    
    //checking all axes for non-intersection of the dimensional interval
    for (int i=0;i<3;i++)
      if ((VMax1(i)<VMin2(i))||(VMax2(i)<VMin1(i)))
        return false;
    
    return true;  //all dimensional intervals are overlapping = intersection
    
  }
  
  bool isCollide(const Mesh& m, double& depth, RowVector3d& intNormal, RowVector3d& intPosition){
    
    
    if ((isFixed && m.isFixed))  //collision does nothing
      return false;
    
    //collision between bounding boxes
    if (!isBoxCollide(m))
      return false;
    
    //otherwise, full test
    ccd_t ccd;
    CCD_INIT(&ccd);
    ccd.support1       = support; // support function for first object
    ccd.support2       = support; // support function for second object
    ccd.center1         =center;
    ccd.center2         =center;
    
    ccd.first_dir       = stub_dir;
    ccd.max_iterations = 100;     // maximal number of iterations
    
    
    void* obj1=(void*)this;
    void* obj2=(void*)&m;
    
    ccd_real_t _depth;
    ccd_vec3_t dir, pos;
    
    int nonintersect = ccdMPRPenetration(obj1, obj2, &ccd, &_depth, &dir, &pos);
    
    if (nonintersect)
      return false;
    
    for (int k=0;k<3;k++){
      intNormal(k)=dir.v[k];
      intPosition(k)=pos.v[k];
    }
    
    depth =_depth;
    intPosition-=depth*intNormal/2.0;
    
    //Vector3d p1=intPosition+depth*intNormal;
    //Vector3d p2=intPosition;
    //std::cout<<"intPosition: "<<intPosition<<std::endl;
    
    //std::cout<<"depth: "<<depth<<std::endl;
    //std::cout<<"After ccdGJKIntersect"<<std::endl;
    
    //return !nonintersect;
    
    return true;
    
  }
  
  
  //return the current inverted inertia tensor around the current COM. Update it by applying the orientation
  Matrix3d getCurrInvInertiaTensor(){
    Matrix3d R=Q2RotMatrix(orientation); // R: Rotation Matrix
    return R.transpose() * invIT * R;
  }
  

  //Update the current position with updated Com
  void updatePositionByCom(){
    if (isFixed)
      return;  //a fixed object is immobile
    for (int i=0;i<currV.rows();i++)
      currV.row(i)<<QRot(origV.row(i), orientation)+COM;
  }

  //Update the current position and orientation by integrating the linear and angular velocities, and update currV accordingly
  //You need to modify this according to its purpose
  void updatePosition(double timeStep){
    //just forward Euler now
    if (isFixed)
      return;  //a fixed object is immobile
    
    COM += timeStep* comVelocity;

    // Rotation
    Eigen::RowVector4d qAngV;
    qAngV << 0, 0.5 * timeStep * angVelocity;
    orientation += QMult(qAngV, orientation);
    orientation.normalize();
    
    for (int i=0;i<currV.rows();i++)
      currV.row(i)<<QRot(origV.row(i), orientation)+COM;
  }
  
  
  //Updating velocity *instantaneously*. i.e., not integration from acceleration, but as a result of a collision impulse from the "impulses" list
  //You need to modify this for that purpose.
  void updateImpulseVelocities(){
    
    if (isFixed){
      comVelocity.setZero();
      currImpulses.clear();
      angVelocity.setZero();
      return;
    }
    
    //update linear and angular velocity according to all impulses
    for (int i=0; i<currImpulses.size(); i++){
      // should use formula on slide 21 of topic 3 to compute the new COM velocity
      Vector3d relativePos = currImpulses[i].first.transpose() - COM.transpose();
      comVelocity += currImpulses[i].second.transpose() / totalMass; // v(t + timeStep) = v(t) + I/m
      angVelocity += getCurrInvInertiaTensor() * relativePos.cross(currImpulses[i].second.transpose()); // w = I^-1 * r x I
    }
    currImpulses.clear();
  }
  
  RowVector3d initStaticProperties(const double density)
  {
    tetVolumes.conservativeResize(T.rows());
    
    RowVector3d naturalCOM; naturalCOM.setZero();
    Matrix3d IT; IT.setZero();
    for (int i=0;i<T.rows();i++){
      Vector3d e01=origV.row(T(i,1))-origV.row(T(i,0));
      Vector3d e02=origV.row(T(i,2))-origV.row(T(i,0));
      Vector3d e03=origV.row(T(i,3))-origV.row(T(i,0));
      Vector3d tetCentroid=(origV.row(T(i,0))+origV.row(T(i,1))+origV.row(T(i,2))+origV.row(T(i,3)))/4.0;
      tetVolumes(i)=std::abs(e01.dot(e02.cross(e03)))/6.0;
      
      naturalCOM+=tetVolumes(i)*tetCentroid;
      
    }
    
    totalVolume=tetVolumes.sum();
    totalMass=density*totalVolume;
    naturalCOM.array()/=totalVolume;
    
    //computing inertia tensor
    for (int i=0;i<T.rows();i++){
      RowVector4d xvec; xvec<<origV(T(i,0),0)-naturalCOM(0),origV(T(i,1),0)-naturalCOM(0),origV(T(i,2),0)-naturalCOM(0),origV(T(i,3),0)-naturalCOM(0);
      RowVector4d yvec; yvec<<origV(T(i,0),1)-naturalCOM(1),origV(T(i,1),1)-naturalCOM(1),origV(T(i,2),1)-naturalCOM(1),origV(T(i,3),1)-naturalCOM(1);
      RowVector4d zvec; zvec<<origV(T(i,0),2)-naturalCOM(2),origV(T(i,1),2)-naturalCOM(2),origV(T(i,2),2)-naturalCOM(2),origV(T(i,3),2)-naturalCOM(2);
      
      double I00, I11, I22, I12, I21, I01, I10, I02, I20;
      Matrix4d sumMat=Matrix4d::Constant(1.0)+Matrix4d::Identity();
      I00 = density*6*tetVolumes(i)*(yvec*sumMat*yvec.transpose()+zvec*sumMat*zvec.transpose()).sum()/120.0;
      I11 = density*6*tetVolumes(i)*(xvec*sumMat*xvec.transpose()+zvec*sumMat*zvec.transpose()).sum()/120.0;
      I22 = density*6*tetVolumes(i)*(xvec*sumMat*xvec.transpose()+yvec*sumMat*yvec.transpose()).sum()/120.0;
      I12 = I21 = -density*6*tetVolumes(i)*(yvec*sumMat*zvec.transpose()).sum()/120.0;
      I10 = I01 = -density*6*tetVolumes(i)*(xvec*sumMat*zvec.transpose()).sum()/120.0;
      I20 = I02 = -density*6*tetVolumes(i)*(xvec*sumMat*yvec.transpose()).sum()/120.0;
      
      Matrix3d currIT; currIT<<I00, I01, I02,
      I10, I11, I12,
      I20, I21, I22;
      
      IT+=currIT;
      
    }
    invIT=IT.inverse();
    if (isFixed)
      invIT.setZero();  //infinite resistance to rotation
  
    return naturalCOM;
    
  }
  
  
  //Updating the linear and angular velocities of the object
  //You need to modify this to integrate from acceleration in the field (basically gravity)
  void updateVelocity(double timeStep, float DragForceCoeff, RowVector3d& OverallAddonComSpeed, RowVector3d& OverallAddonAngSpeed){
    
    if (isFixed)
      return;

    // integrating drag forces
    if (T.rows()>1) 
      dragByTet(DragForceCoeff);
    else
      simpleDrag(DragForceCoeff);

    // integrating gravity
    Vector3d gravity; gravity<<0,-9.8,0.0;
    comVelocity += gravity*timeStep;
    comVelocity += OverallAddonComSpeed;
    angVelocity += OverallAddonAngSpeed;
  }

  // simple drag implemention for object with only 1 tet
  void simpleDrag(float DragForceCoeff){
    comVelocity -= comVelocity * DragForceCoeff / totalMass;
    angVelocity -= angVelocity * DragForceCoeff / totalMass;
  }

  // realistic drag implemention.
  void dragByTet(float DragForceCoeff){
    Vector3d forceCOM; forceCOM << 0,0,0;
    Vector3d torch; torch << 0,0,0; // torque?
    for (int i=0;i<T.rows();i++){
      Vector3d relativePos = getTetCentroid(currV, i) - COM.transpose(); // what is relativePos ? between COM and force contact point?
      precisionCheckZero(relativePos);
      double percent = tetVolumes(i) / totalVolume;
      Vector3d force = -DragForceCoeff * percent * getPointSpeed(relativePos);
      forceCOM += force;
      torch += relativePos.cross(force); // torque = relativePos x force
    }
    precisionCheckZero(forceCOM);
    comVelocity += forceCOM / totalMass; // v1 = v0 + F/m
    angVelocity += getCurrInvInertiaTensor() * torch; // w1 = w0 + I^-1 * torque
  }

  // get the central of the tet
  Vector3d getTetCentroid(const MatrixXd& usedV, int i){
    return (usedV.row(T(i,0))+usedV.row(T(i,1))+usedV.row(T(i,2))+usedV.row(T(i,3)))/4.0;
  }

  // get total speed for a point (user should ensure the point is inside the mesh)
  Vector3d getPointSpeed(const Vector3d& relativePos){
    Vector3d linearAngVelocity = comVelocity;
    Matrix3d angToLinear;
    angToLinear <<  0, -angVelocity(0,2), angVelocity(0,1), 
                    angVelocity(0,2), 0 , -angVelocity(0,0),
                    -angVelocity(0,1), angVelocity(0,0), 0;
    linearAngVelocity += angToLinear * relativePos;
    return linearAngVelocity;
  }

  RowVector3d getPointSpeedRow(const RowVector3d& relativePos){
    return getPointSpeed(relativePos.transpose()).transpose();
  }
  
  
  //the full integration for the time step (velocity + position)
  //You need to modify this if you are changing the integration
  void integrate(double timeStep, float DragForceCoeff, RowVector3d& OverallAddonComSpeed, RowVector3d& OverallAddonAngSpeed){
    updateVelocity(timeStep, DragForceCoeff, OverallAddonComSpeed, OverallAddonAngSpeed);
    updatePosition(timeStep);
  }
  
  void manualAddVelocity(const RowVector3d& AddonComSpeed, const RowVector3d& AddonAngSpeed){
    if (isFixed)
      return;
    comVelocity += AddonComSpeed;
    angVelocity += AddonAngSpeed;
  }
  
  Mesh(const MatrixXd& _V, const MatrixXi& _F, const MatrixXi& _T, const double density, const bool _isFixed, const RowVector3d& _COM, const RowVector4d& _orientation){
    origV=_V;
    F=_F;
    T=_T;
    isFixed=_isFixed;
    COM=_COM;
    orientation=_orientation;
    comVelocity.setZero();
    angVelocity.setZero();
    
    RowVector3d naturalCOM;  //by the geometry of the object
    
    //initializes the original geometric properties (COM + IT) of the object
    naturalCOM = initStaticProperties(density);
    
    origV.rowwise()-=naturalCOM;  //removing the natural COM of the OFF file (natural COM is never used again)
    
    currV.resize(origV.rows(), origV.cols());
    for (int i=0;i<currV.rows();i++)
      currV.row(i)<<QRot(origV.row(i), orientation)+COM;
    
    
    VectorXi boundVMask(origV.rows());
    boundVMask.setZero();
    for (int i=0;i<F.rows();i++)
      for (int j=0;j<3;j++)
        boundVMask(F(i,j))=1;
    
    //cout<<"boundVMask.sum(): "<<boundVMask.sum()<<endl;
    
    vector<int> boundTList;
    for (int i=0;i<T.rows();i++){
      int incidence=0;
      for (int j=0;j<4;j++)
        incidence+=boundVMask(T(i,j));
      if (incidence>2)
        boundTList.push_back(i);
    }
    
    boundTets.resize(boundTList.size());
    for (int i=0;i<boundTets.size();i++)
      boundTets(i)=boundTList[i];
  }
  
  ~Mesh(){}
};

//This class contains the entire scene operations, and the engine time loop.
class Scene{
public:
  double currTime;
  int numFullV, numFullT;
  std::vector<Mesh> meshes;
  float DragForceCoeff;
  float friction;
  
  //Practical 2
  vector<Constraint> constraints;   //The (user) constraints of the scene
  
  //adding an objects. You do not need to update this generally
  void addMesh( const MatrixXd& V, const MatrixXi& F, const MatrixXi& T, const double density, 
                const bool isFixed, const RowVector3d& COM, const RowVector4d& orientation,
                const RowVector3d& InitComVelocity, const RowVector3d& InitAngVelocity){
    
    Mesh m(V,F, T, density, isFixed, COM, orientation);
    m.manualAddVelocity(InitComVelocity, InitAngVelocity);
    meshes.push_back(m);
  }
  
  
  
  /*********************************************************************
   This function handles collision constraints between objects m1 and m2 when found
   Input: meshes m1, m2
   depth: the depth of penetration
   contactNormal: the normal of the conact measured m1->m2
   penPosition: a point on m2 such that if m2 <= m2 + depth*contactNormal, then penPosition+depth*contactNormal is the common contact point
   CRCoeff: the coefficient of restitution
   
   You should create a "Constraint" class, and use its resolveVelocityConstraint() and resolvePositionConstraint() *alone* to resolve the constraint.
   You are not allowed to use practical 1 collision handling
   *********************************************************************/
  void handleCollision(Mesh& m1, Mesh& m2,const double& depth, const RowVector3d& contactNormal,const RowVector3d& penPosition, const double CRCoeff, const double tolerance){
    
    cout <<endl<< "col begin"<<endl;
    
    //std::cout<<"contactNormal: "<<contactNormal<<std::endl;
    //std::cout<<"penPosition: "<<penPosition<<std::endl;
    
    double invMass1 = (m1.isFixed ? 0.0 : 1.0/m1.totalMass);  //fixed meshes have infinite mass
    double invMass2 = (m2.isFixed ? 0.0 : 1.0/m2.totalMass);
  
    /***************
     TODO: practical 2
     update m(1,2) comVelocity, angVelocity and COM variables by using a Constraint class of type COLLISION
     ***********************/

    Constraint collisionConstraint = Constraint(
        COLLISION, INEQUALITY, 0, 1, 0, 1, invMass1, invMass2, contactNormal, 0, CRCoeff);

    Matrix3d invInertiaTensor_m1 = m1.getCurrInvInertiaTensor();
    Matrix3d invInertiaTensor_m2 = m2.getCurrInvInertiaTensor();

    MatrixXd COMPositions(2, 3); COMPositions << m1.COM, m2.COM;
    MatrixXd ConstPositions(2, 3); ConstPositions << penPosition + depth * contactNormal, penPosition;
    MatrixXd COMVelocities(2, 3); COMVelocities << m1.comVelocity, m2.comVelocity;
    MatrixXd AngVelocities(2, 3); AngVelocities << m1.angVelocity, m2.angVelocity;

    MatrixXd correctedCOMVelocities(2, 3), correctedAngVelocities(2, 3), correctedCOMPositions(2, 3);

    bool velocityWasValid = collisionConstraint.resolveVelocityConstraint(
      COMPositions, ConstPositions, COMVelocities, AngVelocities,
      invInertiaTensor_m1, invInertiaTensor_m2,
      correctedCOMVelocities, correctedAngVelocities, tolerance
    );
    bool positionWasValid = collisionConstraint.resolvePositionConstraint(COMPositions, ConstPositions, correctedCOMPositions, tolerance);

    cout << "pos before" <<COMPositions<<endl;
    cout << "pos after" <<correctedCOMPositions<<endl;
    if (velocityWasValid == false) {
        m1.comVelocity = correctedCOMVelocities.row(0);
        m2.comVelocity = correctedCOMVelocities.row(1);
        m1.angVelocity = correctedAngVelocities.row(0);
        m2.angVelocity = correctedAngVelocities.row(1);
    }
    if (positionWasValid == false) {
        m1.COM = correctedCOMPositions.row(0);
        m2.COM = correctedCOMPositions.row(1);
        //updating currV according to corrected COM
        for (int i = 0; i < m1.currV.rows(); i++)
            m1.currV.row(i) << QRot(m1.origV.row(i), m1.orientation) + m1.COM;
        for (int i = 0; i < m2.currV.rows(); i++)
            m2.currV.row(i) << QRot(m2.origV.row(i), m2.orientation) + m2.COM;
    }
  }
  
  /*********************************************************************
   This function handles a single time step by:
   1. Integrating velocities, positions, and orientations by the timeStep
   2. (Practical 2) Detecting collisions and encoding them as constraints
   3. (Practical 2) Iteratively resolved positional and velocity constraints
   
   You do not need to update this function in Practical 2
   *********************************************************************/
  void updateScene( const double timeStep, const double CRCoeff, const double tolerance, 
                    const int maxIterations,float dragForceCoeff, float fricTion, float constraintFlexibleRatio,
                    RowVector3d& OverallAddonComSpeed, RowVector3d& OverallAddonAngSpeed)
{
    DragForceCoeff = dragForceCoeff;
    friction = fricTion;
    for (int i=0; i<constraints.size(); i++)
      constraints[i].ConstraintFlexibleRatio = constraintFlexibleRatio;
    //integrating velocity, position and orientation from forces and previous states
    for (int i=0;i<meshes.size();i++)
      meshes[i].integrate(timeStep, DragForceCoeff, OverallAddonComSpeed, OverallAddonAngSpeed);

    //detecting and handling collisions when found
    //This is done exhaustively: checking every two objects in the scene.
    double depth;
    RowVector3d contactNormal, penPosition;
    for (int i=0;i<meshes.size();i++)
      for (int j=i+1;j<meshes.size();j++)
        if (meshes[i].isCollide(meshes[j],depth, contactNormal, penPosition))
          handleCollision(meshes[i], meshes[j],depth, contactNormal, penPosition,CRCoeff, tolerance);
    
    
    //Resolving user constraints iteratively until either:
    //1. Positions or velocities are valid up to tolerance (a full streak of validity in the iteration)
    //2. maxIterations has run out
    
    
    //Resolving velocity
    int currIteration=0;
    int zeroStreak=0;  //how many consecutive constraints are already below tolerance without any change; the algorithm stops if all are.
    int currConstIndex=0;
    while ((zeroStreak<constraints.size())&&(currIteration*constraints.size()<maxIterations)){
      
      Constraint currConstraint=constraints[currConstIndex];
      
      RowVector3d origConstPos1=meshes[currConstraint.m1].origV.row(currConstraint.v1);
      RowVector3d origConstPos2=meshes[currConstraint.m2].origV.row(currConstraint.v2);
      
      RowVector3d currConstPos1 = QRot(origConstPos1, meshes[currConstraint.m1].orientation)+meshes[currConstraint.m1].COM;
      RowVector3d currConstPos2 = QRot(origConstPos2, meshes[currConstraint.m2].orientation)+meshes[currConstraint.m2].COM;
      //cout<<"(currConstPos1-currConstPos2).norm(): "<<(currConstPos1-currConstPos2).norm()<<endl;
      //cout<<"(meshes[currConstraint.m1].currV.row(currConstraint.v1)-meshes[currConstraint.m2].currV.row(currConstraint.v2)).norm(): "<<(meshes[currConstraint.m1].currV.row(currConstraint.v1)-meshes[currConstraint.m2].currV.row(currConstraint.v2)).norm()<<endl;
      MatrixXd currCOMPositions(2,3); currCOMPositions<<meshes[currConstraint.m1].COM, meshes[currConstraint.m2].COM;
      MatrixXd currConstPositions(2,3); currConstPositions<<currConstPos1, currConstPos2;
      MatrixXd currCOMVelocities(2,3); currCOMVelocities<<meshes[currConstraint.m1].comVelocity, meshes[currConstraint.m2].comVelocity;
      MatrixXd currAngVelocities(2,3); currAngVelocities<<meshes[currConstraint.m1].angVelocity, meshes[currConstraint.m2].angVelocity;
      
      Matrix3d invInertiaTensor1=meshes[currConstraint.m1].getCurrInvInertiaTensor();
      Matrix3d invInertiaTensor2=meshes[currConstraint.m2].getCurrInvInertiaTensor();
      MatrixXd correctedCOMVelocities(2,3), correctedAngVelocities(2,3), correctedCOMPositions(2,3);
      
      bool velocityWasValid=currConstraint.resolveVelocityConstraint(currCOMPositions, currConstPositions, currCOMVelocities, currAngVelocities, invInertiaTensor1, invInertiaTensor2, correctedCOMVelocities,correctedAngVelocities, tolerance);
      
      if (velocityWasValid){
        zeroStreak++;
      }else{
        //only update the COM and angular velocity, don't both updating all currV because it might change again during this loop!
        zeroStreak=0;
        meshes[currConstraint.m1].comVelocity =correctedCOMVelocities.row(0);
        meshes[currConstraint.m2].comVelocity =correctedCOMVelocities.row(1);
        
        meshes[currConstraint.m1].angVelocity =correctedAngVelocities.row(0);
        meshes[currConstraint.m2].angVelocity =correctedAngVelocities.row(1);
        
      }
      
      currIteration++;
      currConstIndex=(currConstIndex+1)%(constraints.size());
    }
    
    if (currIteration*constraints.size()>=maxIterations)
      cout<<"Velocity Constraint resolution reached maxIterations without resolving!"<<endl;
    
    
    //Resolving position
    currIteration=0;
    zeroStreak=0;  //how many consecutive constraints are already below tolerance without any change; the algorithm stops if all are.
    currConstIndex=0;
    while ((zeroStreak<constraints.size())&&(currIteration*constraints.size()<maxIterations)){
      
      Constraint currConstraint=constraints[currConstIndex];
      
      RowVector3d origConstPos1=meshes[currConstraint.m1].origV.row(currConstraint.v1);
      RowVector3d origConstPos2=meshes[currConstraint.m2].origV.row(currConstraint.v2);
      
      RowVector3d currConstPos1 = QRot(origConstPos1, meshes[currConstraint.m1].orientation)+meshes[currConstraint.m1].COM;
      RowVector3d currConstPos2 = QRot(origConstPos2, meshes[currConstraint.m2].orientation)+meshes[currConstraint.m2].COM;

      MatrixXd currCOMPositions(2,3); currCOMPositions<<meshes[currConstraint.m1].COM, meshes[currConstraint.m2].COM;
      MatrixXd currConstPositions(2,3); currConstPositions<<currConstPos1, currConstPos2;
    
      MatrixXd correctedCOMPositions(2,3);
    
      bool positionWasValid=currConstraint.resolvePositionConstraint(currCOMPositions, currConstPositions,correctedCOMPositions, tolerance);
      
      if (positionWasValid){
        zeroStreak++;
      }else{
        //only update the COM and angular velocity, don't both updating all currV because it might change again during this loop!
        zeroStreak=0;

        meshes[currConstraint.m1].COM =correctedCOMPositions.row(0);
        meshes[currConstraint.m2].COM =correctedCOMPositions.row(1);
        
      }
      
      currIteration++;
      currConstIndex=(currConstIndex+1)%(constraints.size());
    }
    
    if (currIteration*constraints.size()>=maxIterations)
      cout<<"Position Constraint resolution reached maxIterations without resolving!"<<endl;
    
    
    //updating currV according to corrected COM
    for (int i=0;i<meshes.size();i++)
      for (int j=0;j<meshes[i].currV.rows();j++)
        meshes[i].currV.row(j)<<QRot(meshes[i].origV.row(j), meshes[i].orientation)+meshes[i].COM;
    
    currTime+=timeStep;
    OverallAddonComSpeed *= 0;
    OverallAddonAngSpeed *= 0;
  }
  
  //loading a scene from the scene .txt files
  //you do not need to update this function
  bool loadScene(const std::string dataFolder, const std::string sceneFileName, const std::string constraintFileName){
    
    ifstream sceneFileHandle, constraintFileHandle;
    sceneFileHandle.open(dataFolder+std::string("/")+sceneFileName);
    if (!sceneFileHandle.is_open())
      return false;
    int numofObjects;
    string line;
    
    currTime=0;
    sceneFileHandle>>numofObjects;
    getline(sceneFileHandle, line);
    for (int i=0;i<numofObjects;i++){
      MatrixXi objT, objF;
      MatrixXd objV;
      std::string MESHFileName;
      bool isFixed;
      double youngModulus, poissonRatio, density;
      RowVector3d userCOM;
      RowVector4d userOrientation;
      sceneFileHandle>>MESHFileName>>density>>youngModulus>>poissonRatio>>isFixed>>userCOM(0)>>userCOM(1)>>userCOM(2)>>userOrientation(0)>>userOrientation(1)>>userOrientation(2)>>userOrientation(3);
      userOrientation.normalize();

      // read optional input, in an order of comv, coma
      RowVector3d InitComVelocity;
      RowVector3d InitAngVelocity;
      getline(sceneFileHandle, line);
      istringstream ss(line);
      if (ss >> InitComVelocity[0])
        ss >> InitComVelocity[1], InitComVelocity[2];
      else
        InitComVelocity.setZero();
      if (ss >> InitAngVelocity[0])
        ss >> InitAngVelocity[1], InitAngVelocity[2];
      else
        InitAngVelocity.setZero();

      igl::readMESH(dataFolder+std::string("/")+MESHFileName,objV,objT, objF);
      
      //fixing weird orientation problem
      MatrixXi tempF(objF.rows(),3);
      tempF<<objF.col(2), objF.col(1), objF.col(0);
      objF=tempF;
      
      addMesh(objV,objF, objT,density, isFixed, userCOM, userOrientation, InitComVelocity, InitAngVelocity);
      cout << "COM: " << userCOM <<endl;
      cout << "orientation: " << userOrientation <<endl;
    }
    
    //Practical 2 change
    //reading intra-mesh attachment constraints
    int numofConstraints;
    constraintFileHandle.open(dataFolder+std::string("/")+constraintFileName);
    if (!constraintFileHandle.is_open())
      return false;
    constraintFileHandle>>numofConstraints;
    for (int i=0;i<numofConstraints;i++){
      int attachM1, attachM2, attachV1, attachV2;
      constraintFileHandle>>attachM1>>attachV1>>attachM2>>attachV2;
      
      double initDist=(meshes[attachM1].currV.row(attachV1)-meshes[attachM2].currV.row(attachV2)).norm();
      //cout<<"initDist: "<<initDist<<endl;
      double invMass1 = (meshes[attachM1].isFixed ? 0.0 : 1.0/meshes[attachM1].totalMass);  //fixed meshes have infinite mass
      double invMass2 = (meshes[attachM2].isFixed ? 0.0 : 1.0/meshes[attachM2].totalMass);
      constraints.push_back(Constraint(DISTANCE, EQUALITY,attachM1, attachV1, attachM2, attachV2, invMass1,invMass2,RowVector3d::Zero(), initDist, 0.0));
      
    }
    
    return true;
  }
  
  
  Scene(){}
  ~Scene(){}
};


/*****************************Auxiliary functions for collision detection. Do not need updating********************************/

/** Support function for libccd*/
void support(const void *_obj, const ccd_vec3_t *_d, ccd_vec3_t *_p)
{
  // assume that obj_t is user-defined structure that holds info about
  // object (in this case box: x, y, z, pos, quat - dimensions of box,
  // position and rotation)
  //std::cout<<"calling support"<<std::endl;
  Mesh *obj = (Mesh *)_obj;
  RowVector3d p;
  RowVector3d d;
  for (int i=0;i<3;i++)
    d(i)=_d->v[i]; //p(i)=_p->v[i];
  
  d.normalize();
  //std::cout<<"d: "<<d<<std::endl;
  
  int maxVertex=-1;
  int maxDotProd=-32767.0;
  for (int i=0;i<obj->currV.rows();i++){
    double currDotProd=d.dot(obj->currV.row(i)-obj->COM);
    if (maxDotProd < currDotProd){
      maxDotProd=currDotProd;
      //std::cout<<"maxDotProd: "<<maxDotProd<<std::endl;
      maxVertex=i;
    }
    
  }
  //std::cout<<"maxVertex: "<<maxVertex<<std::endl;
  
  for (int i=0;i<3;i++)
    _p->v[i]=obj->currV(maxVertex,i);
  
  //std::cout<<"end support"<<std::endl;
}

void stub_dir(const void *obj1, const void *obj2, ccd_vec3_t *dir)
{
  dir->v[0]=1.0;
  dir->v[1]=0.0;
  dir->v[2]=0.0;
}

void center(const void *_obj,ccd_vec3_t *center)
{
  Mesh *obj = (Mesh *)_obj;
  for (int i=0;i<3;i++)
    center->v[i]=obj->COM(i);
}








#endif
