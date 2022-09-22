//
// Created by eecs on 25/02/20.
//

#include <iostream>

using namespace std;

int main (){
  int ma[4][5];
  int a=0;
  for (int i=0; i<4; i++){
    for (int j=0; j<5; j++){
      ma[i][j]=a;
      a++;
    }
    cout << endl;
  }

  for (int i=0; i<4; i++){
    for (int j=0; j<5; j++){
      cout << ma[i][j] << " ";
    }
    cout << endl;
  }


  return 0;
}