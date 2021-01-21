// File:          robot01.cpp
// Date:
// Description:
 // Author:
// Modifications:

// You may need to add webots include files such as
// <webots/DistanceSensor.hpp>, <webots/Motor.hpp>, etc.
// and/or to add some other includes
#include <webots/Robot.hpp>
#include <webots/Motor.hpp>
#include <webots/Keyboard.hpp>
#include <webots/PositionSensor.hpp>
#include <iostream>
#include <cstdint>
#include <cstdio>
#include <math.h>
#include <vector>
#include <initializer_list>

#define MAX_VELOCITY 5
#define MIN_VELOCITY -5
#define MAX_GEAR 5
#define MIN_GEAR -2

// All the webots classes are defined in the "webots" namespace
using namespace webots;
double gear_velocities[8] = {-1,-0.5,0,0.3,0.5,0.7,1,1.2};
double ARMS[3] = {0.45, 0.1, 0.5};
const double PI = 3.141592;
enum Side {RIGHT,LEFT,NONE};


class Matrix{
private:
    unsigned int rows;
    unsigned int columns;
    std::vector<std::vector<double> > values;
    
   double mac(std::vector<double> row, std::vector<double> column){
    	double sum = 0;
		for(unsigned int i = 0; i < this->rows ; i++ ){
			sum += row.at(i) * column.at(i);
			std::cout << row.at(i) << "\t" << column.at(i) << std::endl;
		}
		return sum;
	}
	//END OF MATRIX::MULTIPLY-AND-ACCUMULATE

public:

    Matrix(){
      rows = 3;
      columns = 3;
      values.resize(3);
        for(unsigned int r = 0; r < 3; r++){
            values.at(r).resize(3);
        }
    }

    Matrix(unsigned int rows, unsigned int columns){
        if(rows < 1)
            rows = 1;
        if(columns < 1)
            columns = 1;

        this->columns = columns;
        this->rows = rows;
        values.resize(rows);
        for(unsigned int r = 0; r < rows; r++){
            values.at(r).resize(columns);
        }
    }
    //END OF MATRIX::CONSTRUCTOR

    Matrix operator*(Matrix mat){
		if(this->columns != mat.rows)
			throw;
		Matrix res(this->rows, mat.columns);
		
		for(unsigned int r = 0; r < res.rows ; r++){
			for(unsigned int c = 0; c < res.columns ; c++){
				res.get(r,c) = res.mac(values.at(r),mat.get_column(c));
			}	
		}
		return res;
    }
	//END MATRIX::OPERATOR*
	
  std::vector<double> get_column(unsigned int column){
    std::vector<double> res;
    res.resize(this->columns);
    for(unsigned int r = 0 ; r < this->rows ; r++){
      res.at(r) = this->values.at(r).at(column);
    }	
    return res;
    }
    //END MATRIX::GET_COLUMN
    Matrix operator!(){
        Matrix transpose(columns,rows);
        for(int r = 0 ; r < rows ; r++){
            for(int c = 0 ; c < rows ; c++){
                transpose.get(c,r) = values.at(r).at(c);
            }
        }
        return transpose;
    }
    //END MATRIX::OPERATOR! (TRANSPOSE)
    Matrix operator^(int n){

    }
    //END MATRIX::OPERATOR^

    double& get(unsigned int row, unsigned int col){
        return (values.at(row).at(col));
    }
    //END MATRIX::OPERATOR[]

    void print_values(){
        for(int r = 0 ; r < rows; r++){
            for(int c = 0 ; c < columns ;  c++){
                std::cout << values.at(r).at(c) << "  ";
            }
            std::cout << std::endl;
        }
    }
    
    void set(std::vector<double> vals){
		unsigned int counter = 0;
		for(unsigned int r = 0; r < rows ; r++){
			for(unsigned int c = 0; c < columns ; c++){
				values.at(r).at(c) = vals[counter];
				counter++;
			}
		}
	}

};

// GRIPPER CLASS
class Gripper{
  public:
    Gripper(Motor **fingers){ // right finger first
      this->fingers = fingers;
    }
    
    void set_velocities(){
      for(int i = 0; i < 2; i++){
        fingers[i]->setVelocity(0.1);
      }
    }
    
    bool is_open(){}
    bool is_closed(){}
    
    void close(){
      fingers[0]->setPosition(fingers[0]->getMaxPosition());
      fingers[1]->setPosition(fingers[1]->getMinPosition());
     this->set_velocities();
    }

    void open(){
       fingers[0]->setPosition(fingers[0]->getMinPosition());
       fingers[1]->setPosition(fingers[1]->getMaxPosition());
       this->set_velocities();
    }
    
    void release(){}
    
    void tighten(){
      double right_finger_pos = fingers[0]->getPositionSensor()->getValue();
      double left_finger_pos = fingers[1]->getPositionSensor()->getValue();
      std::cout << right_finger_pos << std::endl << left_finger_pos<<std::endl;
      fingers[0]->setPosition(right_finger_pos+0.01);
      fingers[1]->setPosition(left_finger_pos-0.01);
      this->set_velocities();
    }

  private:
    Motor **fingers;
};

//ENGINE CLASS
class Engine{
  public:
    Engine(Motor **tracks){ // right track first
      this->tracks = tracks;
      this->gear = 0;
    }
    
    void stop(){
      gear = 0;
      rotation = NONE;
    }
    
    void rotate(Side side){
      this->rotation = side;
    }
    void change_gear(bool up){
      if(up && this->gear < 5) 
        gear++;
      else if(!up && this->gear > -2) 
        gear--;
    }
    void turn(Side side){
      turn_side = side;
    }
    int get_gear(){}
    void set_gear(){}
    
    bool is_running(){
      if (gear != 0){
        return true;
      }
      return false;
    return false;
    }
    
    void run(){
      if(rotation == RIGHT){
        tracks[0]->setVelocity(gear_velocities[gear+2]);
        tracks[1]->setVelocity((-1)*gear_velocities[gear+2]);
      }
      else if(rotation == LEFT){
        tracks[0]->setVelocity((-1)*gear_velocities[gear+2]);
        tracks[1]->setVelocity(gear_velocities[gear+2]);
      }
      else if(turn_side == RIGHT){
        tracks[0]->setVelocity(gear_velocities[gear+2]);
        tracks[1]->setVelocity(gear_velocities[gear+2]/2);
      }
      else if(turn_side == LEFT){
        tracks[0]->setVelocity(gear_velocities[gear+2]/2);
        tracks[1]->setVelocity(gear_velocities[gear+2]);
      }
      else{
        tracks[0]->setVelocity(gear_velocities[gear+2]);
        tracks[1]->setVelocity(gear_velocities[gear+2]);
      }
      
    }
  private:
  Motor **tracks;
  int gear = 0;
  Side rotation = NONE;
  Side turn_side = NONE;
};

class Arm{
  public:
    Arm(Motor **joints){
      this->joints = joints;
    }
    //END OF CONSTRUCTOR
    //BEGIN RESET FUNCTION
    void reset(){
      for(int i = 0; i < 4; i++){
        if(i==0)
          joints[i]->setPosition(recalculate_base_position());
        else
          joints[i]->setPosition(reset_positions[i]);
        joints[i]->setVelocity(1);
        std::cout << joints[i]->getPositionSensor()->getValue() << std::endl;
      }
    }
    //END OF RESET FUNCTION
    //BEGIN ROTATE FUNCTION
    void rotate(Side side){
      joints[0]->setPosition(INFINITY);
      if(side==RIGHT)
        joints[0]->setVelocity(-1.0);
      else if(side == LEFT)
        joints[0]->setVelocity(1.0);
      else
        joints[0]->setVelocity(0.0);
    }
    //END ROTATE FUNCTION
    
    void stop(){
      for(int i = 0; i < 4; i++){
        joints[i]->setVelocity(0.0);
      }
    }
    
    void move_depth(int i){
    return;
      if(i==0)
        return; 
      double qn[3];
      print_end_effector();
      inverse_kinematics(true, qn);
      for(int i = 0; i < 3; i++){
        joints[i+1]->setPosition(qn[i]);
        joints[i+1]->setVelocity(1);
       // std::cout << "Target " << i << "\t" << joints[i+1]->getTargetPosition() << "\t" << qn[i] << std::endl;
       // std::cout << "Current " << i << "\t" << joints[i+1]->getPositionSensor()->getValue() <<std::endl;
        
      }
      print_end_effector();
      }

    void move_height(bool heigher){}
    
    bool control_position(unsigned int index){
      if(joints[index]->getPositionSensor()->getValue() > depth_positions[index])
        return true;
      return false;
    }
    
    double recalculate_base_position(){
      int i = 0;
      double factor = 0;
      double base_value = joints[0]->getPositionSensor()->getValue();
      if (base_value >= 0){
        while(true){
          factor = i*PI;
          if (factor - base_value < 0)
            break;
          i++;
        }
        return 2*PI*i-base_value;
      }
      else{
        while(true){
          factor = i*PI;
          if (factor + base_value > 0)
            break;
          i--;
        }
        return 2*PI*i*(-1)+base_value;
      }
    }
    
    void forward_kinematics(){
      double* q = get_angles();
       end_effector[0] = sin(q[0])*ARMS[0] + ssin(q,0,1)*ARMS[1] + sin(q[2])*ARMS[2];
       end_effector[1] = cos(q[0])*ARMS[0] - scos(q,0,1)*ARMS[1] + +cos(q[2])*ARMS[2];
      delete[] q;
    }
    
    void inverse_kinematics(bool is_depth,double* qn){
      //STEP 1 FIND P2(ie the end of second arm and the begining of the last one)
      double p2[2];
      double yn,xn = 0;
      double theta,d = 0;
      double q1n,q2n,q3n = 0;
      double* q = get_angles();
      double gamma = q[0]+q[1]+q[2];
      get_p(2,p2);
      if(is_depth){
        xn = end_effector[0] + 0.1;
        yn = end_effector[1];
      }
      else{
        xn = end_effector[0];
        yn = end_effector[1] + 0.1;
      }
      d = sqrt(xn*xn+yn*yn);
      theta = atan2((-yn)/d,(-xn)/d);
      q1n = theta + acos((-1)*(xn*xn+yn*yn*ARMS[0]*ARMS[0]-ARMS[1]*ARMS[1])/(d*2*ARMS[0]));
      print_angle(theta,"Theta");
      q2n = atan2((xn-ARMS[0]*cos(q1n))/ARMS[1] , (yn-ARMS[0]*sin(q1n))/ARMS[1] ) - q1n;
      q3n = gamma - q1n - q2n;
      qn[0] = q1n;
      qn[1] = q2n;
      qn[2] - q3n;
      print_angle(q1n,"Q1N");
      print_angle(q[0],"Q1");
      print_angle(q2n,"Q2N");
      print_angle(q[1],"Q2");
      print_angle(q3n,"Q3N");
      print_angle(q[2],"Q3");
      //STEP 2 SOLVE FOR RR PLANAR
      //STEP 3 ???
      //STEP 4 PROFIT
      delete[] q;
    
    }
    
    double rad_to_deg(double rad){
      return rad*180/PI;
    }
    
    void print_angle(double rad,std::string name){
      std::cout << name << "\t(" << rad_to_deg(rad) << " DEG, " << rad << "RAD" << std::endl; 
    }
    
    
    void forward_kinematics_arm(int joint){
    double* q = get_angles();
    Matrix Rs[3];
    Matrix Ts[3];
    int  k = 1;
    for(int i = 0; i < 3; i++){
      if(i==2)
        k = -1;
      Rs[i].set({ -k*cos( sangle(q,0,i ) ), sin( sangle(q,0,i) ), 0,
              sin( sangle(q,0,i) ), -k*cos( sangle(q,0,i) ), 0,
              0, 0, 1 });
      Ts[i].set({1,0,ARMS[i],
            0,1,0,
            0,0,1});
    }
    std::cout << "Gamma = " << (q[0]+q[1]+q[2])*180/PI << " DEG" << std::endl;
     if(joint == 1 || joint == 4){
      std::cout << "ARM1( x= " << sin(q[0])*ARMS[0] << " , y=" << cos(q[0])*ARMS[0] << std::endl;
      (Rs[0]*Ts[0]).print_values();
      }
    else if(joint == 2 || joint == 4 ){
       std::cout << "ARM2( x= " << sin(q[0])*ARMS[0]+ssin(q,0,1)*ARMS[1] << " , y=" << cos(q[0])*ARMS[0]-scos(q,0,1)*ARMS[1] << std::endl;
       (Rs[0]*Ts[0]*Rs[1]*Ts[1]).print_values();
     }
    //else if(joint == 2 || joint == 4 )
       //std::cout << "ARM2( x= " << sin(q[1])*ARMS[1] << " , y=" << cos(q[1])*ARMS[1] << std::endl;
    else if(joint == 3 || joint == 4){
       std::cout << "ARM3-Abs( x= " << sin(q[0])*ARMS[0]+ssin(q,0,1)*ARMS[1]+sin(q[2])*ARMS[2] << " , y=" << cos(q[0])*ARMS[0]-scos(q,0,1)*ARMS[1] +cos(q[2])*ARMS[2] << std::endl;
       std::cout << "ARM3-Rel( x= " << sin(q[2])*ARMS[2] << " , y=" << cos(q[2])*ARMS[2] << std::endl;
       (Rs[0]*Ts[0]*Rs[1]*Ts[1]*Rs[2]*Ts[2]).print_values();
       }
     delete[] q;

    }
        
    void rotate_joint(unsigned int joint_number, Side side){
          int multi = 0;
      if (joint_number > 3 || joint_number < 0){
        return;
      }
      if(side == NONE){
        joints[joint_number]->setVelocity(0);
        return;
      }
      else if(side == RIGHT)
        multi = 1;
      else
        multi = -1;
      forward_kinematics_arm(joint_number);
      //print_end_effector();
      joints[joint_number]->setPosition(INFINITY);
      joints[joint_number]->setVelocity(multi);
    }
    
    double scos(double qs[], int begin, int end){
    // COSINE OF SUMS
      double res = 0;
      for(int i = begin; i < end+1 ; i++)
        res += qs[i];
      return cos(res);
    }
    
    double ssin(double qs[], int begin, int end){
    // COSINE OF SUMS
      double res = 0;
      for(int i = begin; i < end+1 ; i++)
        res += qs[i];
      return sin(res);
    }
    
    double sangle(double qs[], int begin, int end){
      double res = 0;
      for(int i = begin; i < end+1 ; i++)
        res += qs[i];
      return res;
    }
    
    void print_end_effector(){
      std::cout << "x=" << end_effector[0] << " y=" << end_effector[1] << std::endl;
    }
    
    double* get_angles(){
      double* q = new double[3];
      for( int i = 0; i < 3; i++ ){
        q[i] = joints[i+1]->getPositionSensor()->getValue();
        if(i==0)
        q[i] += PI/2-0.26;
        if(i==1){
          q[i] *= (-1);
          q[i] += 0.854;
             }
          if(i==2){
          //q[i] *= (-1);
          q[i] += PI*3/4;
          } 
       }
       return q;
      }
      
     void get_p(unsigned int n, double* result){
       double* q = get_angles();
       if(n==1){
         result[0] = sin(q[0])*ARMS[0];
         result[1] = cos(q[0])*ARMS[0];
         return;
       }
       else if(n==2){
         result[0] = sin(q[0])*ARMS[0]+ssin(q,0,1)*ARMS[1];
         result[1] = cos(q[0])*ARMS[0]-scos(q,0,1)*ARMS[1];
       }
     }
    
  private:
    Motor **joints;;
    const double reset_positions[4] = {0.0, -0.5, 0.0, -0.75};
    double max_positions[4] = {INFINITY, 1.72, 0.5, 0.33};
    double min_positions[4] = {INFINITY, -1.4, -4, -2.73};
    double depth_positions[4] = {0.0, 0.30, -0.5, 0.33};
    double end_effector[2] = {0.0, 0.0}; // depth and height

};

 double max(double a,double b){
   if(a > b) return a;
   return b;
 }
 
 double min(double a, double b){
   if(a < b) return a;
   return b;
 }

 
 void init_device(Robot* robot,Motor** array, bool sensors, int size,std::string *names){
   for (int i = 0; i < size; i++){
     array[i] = robot->getMotor(names[i]);
     array[i]->setPosition(INFINITY);
     array[i]->setVelocity(0.0);
     if(sensors)
       array[i]->getPositionSensor()->enable(int(robot->getBasicTimeStep()/4));
   }
 }
  

// This is the main program of your controller.
// It creates an instance of your Robot instance, launches its
// function(s) and destroys it at the end of the execution.
// Note that only one instance of Robot should be created in
// a controller program.
// The arguments of the main function can be specified by the
// "controllerArgs" field of the Robot node
int main(int argc, char **argv) {
  // CREATE INSTANCES OF WEBOTS NODES/OBJECTS
  Keyboard kb;
  Robot *robot = new Robot();
  int timeStep = (int)robot->getBasicTimeStep();
  Motor *wheels[2];
  Motor *fingers[2];
  Motor *joints[4];
  
  
  std::string fingers_names[2] = {"right_finger","left_finger"};
  std::string wheels_names[2] = {"right","left"}; // "left/0"
  std::string joints_names[4] = {"base","joint1","joint2","joint3"};
    

  init_device(robot,wheels,true,2,wheels_names);
  init_device(robot,fingers,true,2,fingers_names);
  init_device(robot,joints,true,4,joints_names);
  
  Gripper grip(fingers);
  Engine engine(wheels);
  Arm arm(joints);
  std::cout << "XD";
    
  // get the time step of the current world.

  kb.enable(timeStep);
  // You should insert a getDevice-like function in order to get the
  // instance of a device of the robot. Something like:
  //  Motor *motor = robot->getMotor("motorname");
  //  DistanceSensor *ds = robot->getDistanceSensor("dsname");
  //  ds->enable(timeStep);

  // Main loop:
  // - perform simulation steps until Webots is stopping the controller
  while (robot->step(timeStep) != -1) {
    
    int key = kb.getKey();
    if(key == kb.UP)
      engine.change_gear(true);
    else if(key == kb.DOWN)
       engine.change_gear(false);
      
    if(key == kb.RIGHT)
      engine.turn(RIGHT);
    else if(key == kb.LEFT)
      engine.turn(LEFT);
    else
      engine.turn(NONE);
      
     if(key == 'A')
     engine.rotate(LEFT);
     else if(key == 'D')
     engine.rotate(RIGHT);
     else
     engine.rotate(NONE);
    
    
    if(key == 'P')// press p to stop
      engine.stop();
    
    if(key == 'W'){
      arm.move_depth(1);
      std::cout << joints[1]->getTargetPosition() << std::endl;
      }
    else if(key == 'S'){
      arm.move_depth(2);
    }
    else
      arm.move_depth(0);
      
     if(key == 'R')
      arm.reset();
      
         if(key == 'E'){
      arm.rotate(RIGHT);
      }
    else if(key == 'Q')
      arm.rotate(LEFT);
    else
      arm.rotate(NONE);
      
     if(key == '1')
       arm.rotate_joint(1,RIGHT);
     else if(key == '2')
       arm.rotate_joint(1,LEFT);
     else
       arm.rotate_joint(1,NONE);
       
       
       
      if(key == '3')
       arm.rotate_joint(2,RIGHT);
     else if(key == '4')
       arm.rotate_joint(2,LEFT);
     else
       arm.rotate_joint(2,NONE);
       
       
       
      if(key == '5')
       arm.rotate_joint(3,RIGHT);
     else if(key == '6')
       arm.rotate_joint(3,LEFT);
     else
       arm.rotate_joint(3,NONE);
     
      

    if(key == 'G'){
       grip.close();
     }
     else if(key == 'O'){
       grip.open();
     }
     else if(key == 'T'){
       grip.tighten();
     }
     


      
    // Enter here functions to send actuator commands, like:
    //  motor->setPosition(10.0);
    //std::cout << "Q1=" << joints[1]->getPositionSensor()->getValue() << std::endl;1
    arm.forward_kinematics();
    engine.run();
    
  }

  // Enter here exit cleanup code.

  delete robot;
  delete[] joints;
  delete[] wheels;
  delete[] fingers;
  return 0;
}
