#include <iostream>

int main() {
    boolfunc a = boolfunc(vector<bool>{0,0,0,1,1,1,1,0});
    a.ZhegalkinPolynomial();
    cout << "the Zhegalkin Polynomial: " << a.getExp() << endl;
    a.DNF();
    cout << "Perfect DNF: " << a.getExp() << endl;
    a.KNF();
    cout << "Perfect KNF: " << a.getExp() << endl;
    a.BlakeAlg(1);
    cout << "Blake's algorithm (DNF): " << a.getExp() << endl;
    a.BlakeAlg(0);
    cout << "Blake's algorithm (KNF): " << a.getExp() << endl;
    return 0;
}
