#include <iostream>
#include "boolfunc.h"

int main(){
    boolfunc a = boolfunc(vector<bool>{1,1,1,1,0,0,1,0});
    a.BlakeAlg(0);
    return 0;
}
