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

#define MAX_VELOCITY 5
#define MIN_VELOCITY -5
#define MAX_GEAR 5
#define MIN_GEAR -2

// All the webots classes are defined in the "webots" namespace
using namespace webots;
double gear_velocities[8] = {-1,-0.5,0,0.3,0.5,0.7,1,1.2};
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
  enum Side {RIGHT,LEFT,NONE};
  
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
      int multiplicator = 0;
      double left_track = tracks[1]->getPositionSensor()->getValue();
      double right_track = tracks[0]->getPositionSensor()->getValue();
      if(gear < 0){
        rotation = NONE;
  
        return;
      }
      else if(rotation == RIGHT)
        multiplicator = 1;
      else
        multiplicator = -1;
      tracks[1]->setPosition(left_track-0.05*gear*multiplicator);
      tracks[0]->setPosition(right_track+0.05*gear*multiplicator);
    }
    void change_gear(bool up){
      if(up && this->gear < 5) 
        gear++;
      else if(!up && this->gear > -2) 
        gear--;
    }
    void turn(Side side){}
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
      if(tracks[0]->getPositionSensor()->getValue()+0.01 >= tracks[0]->getTargetPosition() && rotation == RIGHT)
        rotation = NONE;
      if (rotation == NONE || gear < 0){
        tracks[1]->setPosition(INFINITY);
        tracks[0]->setPosition(INFINITY);
      }
      
      tracks[0]->setVelocity(gear_velocities[gear+2]);
      tracks[1]->setVelocity(gear_velocities[gear+2]);
      
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
    void reset(){
    
    }
    void set_endeffector_position(){
    
    }
    void rotate(){
      
    }
    
    void move_y_axis(){
    
    }
    void move_x_axis(){
    
    }
  private:
    Motor **joints;

};

 double max(double a,double b){
   if(a > b) return a;
   return b;
 };
 
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
  double velocity_left = 0.0;
  double velocity_right = 0.0;
  double velocity_arm1 = 0.0;
  double velocity_joint1 =0.0;
  double velocity_joint2 =0.0;
  double velocity_joint3 =0.0;
  unsigned int turn = 0;
  unsigned int arm1_turn = 0;
  unsigned int joint1_turn = 0;
  unsigned int joint2_turn = 0;
  unsigned int joint3_turn = 0;
  double turn_time_start = 0;
  double allowed_turn_time = 0.010;
  double arm1_time_start = 0;
  double joint1_time_start = 0;
  double joint2_time_start = 0;
  double joint3_time_start = 0;
  int gear = 0;
    

  init_device(robot,wheels,true,2,wheels_names);
  init_device(robot,fingers,true,2,fingers_names);
  init_device(robot,joints,false,4,joints_names);
  
  Gripper grip(fingers);
  Engine engine(wheels);
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
    if(robot->getTime() >= turn_time_start + allowed_turn_time && turn != 0)
      turn = 0;
    
    if(robot->getTime() >= arm1_time_start + allowed_turn_time && arm1_turn != 0)
      arm1_turn = 0;  
    
    if(robot->getTime() >= joint1_time_start + allowed_turn_time && joint1_turn != 0)
      joint1_turn = 0;
      
     if(robot->getTime() >= joint2_time_start + allowed_turn_time && joint2_turn != 0)
      joint2_turn = 0;
      
     if(robot->getTime() >= joint3_time_start + allowed_turn_time && joint3_turn != 0)
      joint3_turn = 0;
     
          
    int key = kb.getKey();
    if(key == kb.UP){
      gear++;
      if(gear > MAX_GEAR)
        gear = MAX_GEAR;
      std::cout << "XD" << std::endl;
      engine.change_gear(true);
    }
    else if(key == kb.RIGHT){
      engine.rotate(engine.RIGHT);
      turn_time_start = robot->getTime();
    }
    else if(key == kb.LEFT){
      engine.rotate(engine.LEFT);
      turn_time_start = robot->getTime();
    }
    else if(key == kb.DOWN){
       gear--;
       if(gear < MIN_GEAR)
          gear = MIN_GEAR;
       engine.change_gear(false);
    }
    else if(key == 'P'){ // press s to stop
      engine.stop();
    }
    
    else if(key == 'E'){
      arm1_turn = 2;
      arm1_time_start = robot->getTime();
    }
    
    else if(key == 'Q'){
      arm1_turn = 1;
      arm1_time_start = robot->getTime();
    }
    
     else if(key == 'S'){
      joint1_turn = 2;
      joint1_time_start = robot->getTime();
    }
    
    else if(key == 'W'){
      joint1_turn = 1;
      joint1_time_start = robot->getTime();
      }
      
        else if(key == '1'){
      joint2_turn = 2;
      joint2_time_start = robot->getTime();
    }
    
    else if(key == '2'){
      joint2_turn = 1;
      joint2_time_start = robot->getTime();
      }
         else if(key == '3'){
      joint3_turn = 2;
      joint3_time_start = robot->getTime();
    }
    
    else if(key == '4'){
      joint3_turn = 1;
      joint3_time_start = robot->getTime();
      }
    else if(key == 'G'){
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
    engine.run();
    
  };

  // Enter here exit cleanup code.

  delete robot;
  delete[] joints;
  delete[] wheels;
  delete[] fingers;
  return 0;
}
