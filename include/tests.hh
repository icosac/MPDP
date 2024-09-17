/**
* @file tests.hh
* @author Enrico Saccon <enricosaccon96@gmail.com>
* @license This project is released under the GNU Public License 3.0.
* @copyright Copyright 2020 Enrico Saccon. All rights reserved.
* @brief File containing tests run for (Frego et al., 2020).
*/

#pragma once

// System includes
#include<fstream>

// Library includes
#ifndef CUDA_ON
#include<utils.hh>
#else
#include<utils.cuh>
#endif

class Run{
public:
  std::string name, test_name, file_name, guessInitialAngles;
  double time, length, err, initTime, endTime;
  uint discr, threads, functionType, jump, rip;

  Run ( std::string _name, uint _discr, double _time, double _length, double _err,
        std::string _test_name, uint _rip, uint _threads, uint _functionType, 
        uint _jump, std::string  _guessInitialAngles, std::string _file_name="", 
        double _initTime=0, double _endTime=0) 
      : name(_name), discr(_discr), time(_time), length(_length), err(_err), test_name(_test_name), rip(_rip), 
        threads(_threads), functionType(_functionType), jump(_jump), guessInitialAngles(_guessInitialAngles), 
        file_name(_file_name), initTime(_initTime), endTime(_endTime)
  {}

  void write(std::fstream& out){
    out << 
      "{\"name\" : \"" << name << 
      "\", \"test_name\" : \"" << test_name <<  "\"" <<
      ", \"discr\" : " << discr << 
      ", \"time\" : " << time << 
      ", \"length\" : " << std::setprecision(20) << length << 
      ", \"err\" : " << std::setprecision(1) << std::scientific << err << 
      ", \"refinements\" : " << rip <<
      ", \"threads\" : " << threads <<
      ", \"functionType\" : " << functionType <<
      ", \"jump\" : " << jump <<
      ", \"guessInitialAngles\" : \"" << guessInitialAngles <<  "\"" << 
      (file_name=="" ? "" : ", \"power_file\" : \""+file_name+"\"") << 
      ", \"initTime\" : " << std::setprecision(20) << initTime << 
      ", \"endTime\" : " << std::setprecision(20) << endTime << 
      "},\n";
  }
};

