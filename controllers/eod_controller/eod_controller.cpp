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
const double PI = 3.141592;
enum Side {RIGHT,LEFT,NONE};
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
    void set_endeffector_position(){
    
    }
    void rotate(Side side){
      joints[0]->setPosition(INFINITY);
      if(side==RIGHT)
        joints[0]->setVelocity(-1.0);
      else if(side == LEFT)
        joints[0]->setVelocity(1.0);
      else
        joints[0]->setVelocity(0.0);
    }
    
    void stop(){
      for(int i = 0; i < 4; i++){
        joints[i]->setVelocity(0.0);
      }
    }
    
    void move_y_axis(){
      
    
    }
    void move_x_axis(){
    
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
    
  private:
    Motor **joints;;
    const double reset_positions[4] = {0.0, -0.5, 0.0, -0.75};
    const double max_positions[4] = {INFINITY, 1.72, 0.5, 0.33};
    const double min_positions[4] = {0.0, -1.4, -4, -2.73};

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
    
    
    if(key == 'P')// press s to stop
      engine.stop();
    
    if(key == 'E')
      arm.rotate(RIGHT);
    else if(key == 'Q')
      arm.rotate(LEFT);
    else
      arm.rotate(NONE);
      
    if(key == 'R')
      arm.reset();
    
     if(key == 'S'){
    }
    
    else if(key == 'W'){

      }
      
    else if(key == '1'){
    }
    
    else if(key == '2'){
      }
    else if(key == '3'){
    }
    
    else if(key == '4'){
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
