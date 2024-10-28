#include "boolfunc.h"

boolfunc::boolfunc()
{
    exp = "NONE";
    K0_D1 = -1;
    inexp = vector<vector<int>>(2,vector<int>(1,0));
    func = vector<bool>(2,0);
    varnames.clear();
    varnames.push_back("x1");
    idF = 1;
    idT = 0;
    countvars = 1;
}

boolfunc::boolfunc(int vars, bool id_T_F)
{
    func = vector<bool>(1<<vars,id_T_F);
    countvars = vars;

    varnames = vector<string>(vars,"x");
    for(int i = 1; i <= vars; i++)
        varnames[i-1]+=to_string(i);

    exp = "NONE";
    inexp = vector<vector<int>>(1<<vars,vector<int>(vars,0));

    idF = !id_T_F;
    idT = id_T_F;

    K0_D1 = -1;
}

boolfunc::boolfunc(vector<bool> f)
{
    int s = f.size();

    exp = "NONE";
    K0_D1 = -1;

    if(!s or s == 1){
        boolfunc();
        return;
    }

    int counter = 1;

    while(1<<counter < s)
        counter++;

    s = 1<<counter;
    countvars = counter;

    func.clear();
    idT = 1;
    idF = 0;
    for(int i = 0; i < s; i++){
        func.push_back(f[i]);
        idT *= f[i];
        idF += f[i];
    }
    idF = !idF;

    varnames = vector<string>(counter,"x");
    for(int i = 1; i <= counter; i++)
        varnames[i-1]+=to_string(i);

    inexp = vector<vector<int>>(s,vector<int>(counter,0));
}

bool boolfunc::getBit(int val, int n){
    bool q;
    n++;
    while(n){
        q = (((val>>1)<<1)!=val);
        val = val>>1;
        n--;
    }
    return q;
}

void boolfunc::KNF()
{
    if(idT){
        exp = "NONE";
        inexp = vector<vector<int>>(1<<countvars,vector<int>(countvars,0));
        return;
    }
    K0_D1 = 0;
    string e = "";

    int size = 1<<countvars;

    for(int i = 0; i < size; i++){
        if(!func[i]){
            string ine = "(";
            for(int j = 0; j < countvars; j++){
                inexp[i][j] = j+1;
                if(getBit(i,countvars-j-1)){
                    ine+='!';
                    inexp[i][j] *= -1;
                }
                ine+=varnames[j];
                if(j+1 != countvars)
                    ine+=" v ";
            }
            ine+=')';
            e+=ine;
            e+=" ^ ";
        } else {
            inexp[i] = vector<int>(countvars,0);
        }
    }

    if(e[e.size()-2] == '^'){
        e.pop_back();
        e.pop_back();
    }

    exp = e;
}

void boolfunc::DNF()
{
    if(idF){
        exp = "NONE";
        inexp = vector<vector<int>>(1<<countvars,vector<int>(countvars,0));
        return;
    }

    K0_D1 = 1;
    string e = "";

    int size = 1<<countvars;

    for(int i = 0; i < size; i++){
        if(func[i]){
            string ine = "(";
            for(int j = 0; j < countvars; j++){
                inexp[i][j] = j+1;
                if(!getBit(i,countvars-j-1)){
                    ine+='!';
                    inexp[i][j] *= -1;
                }
                ine+=varnames[j];
                if(j+1 != countvars)
                    ine+=" ^ ";
            }
            ine+=')';
            e+=ine;
            e+=" v ";
        } else {
            inexp[i] = vector<int>(countvars,0);
        }
    }

    if(e[e.size()-2] == 'v'){
        e.pop_back();
        e.pop_back();
    }

    exp = e;
}

bool boolfunc::Resolution()
{
    int x = inexp.size();

    for(int i = 0; i < x; i++){
        for(int j = 0; j < x; j++){
            int c = 0;
            for(int val = 0; val < countvars; val++){
                if(!(inexp[i][val]+inexp[j][val]) and inexp[i][val] and i != j)
                    c++;
            }

            if(c == 1){
                vector<int> R;

                for(int val = 0; val < countvars; val++)
                    R.push_back(!inexp[i][val] or !inexp[j][val]?!inexp[i][val]?inexp[j][val]:inexp[i][val]:(inexp[i][val]+inexp[j][val])/2);

                bool flag = 1;

                for(int pr = 0; pr < x; pr++)
                    flag *= R != inexp[pr];

                if(flag){
                    string e = "";

                    if(K0_D1)
                        e+=" v (";
                    else
                        e+=" ^ (";

                    int sum = 0;

                    for(int val = 0; val < countvars; val++)
                        sum+=(R[val]!=0);

                    sum--;

                    for(int val = 0; val < countvars; val++){
                        if(R[val]){
                            if(R[val]!=abs(R[val]))
                                e+='!';
                            e+=varnames[val];
                            if(sum){
                                if(K0_D1)
                                    e+=" ^ ";
                                else
                                    e+=" v ";
                                sum--;
                            }
                        }
                    }

                    e+=")";

                    exp+=e;
                    inexp.push_back(R);
                    return 1;
                }
            }
        }
    }
    return 0;
}

void boolfunc::MaxAbsorption()
{   
    vector<int> zero(countvars,0);
    vector<bool> temp(inexp.size(),0);
    reverse(inexp.begin(),inexp.end());

    int x = inexp.size();

    for(int i = 0; i < x-1; i++){
        for(int j = i+1; j < x; j++){
            int sum = 0;
            int countnotnull = 0;
            for(int val = 0; val < countvars; val++)
            {
                countnotnull+=(inexp[i][val])?1:0;
                sum+=(inexp[i][val])?(abs(inexp[i][val]-inexp[j][val])):0;
            }
            temp[j] = temp[j]+(!sum and countnotnull);
        }
    }

    for(int i = 0; i < x; i++)
        if(temp[i])
            inexp[i] = zero;

    x = inexp.size();

    exp = "";
    for(int i = 0; i < x; i++){
        string e = "(";
        int count = 0;
        for(int j = 0; j < countvars; j++)
            count+=!!inexp[i][j];
        for(int j = 0; j < countvars; j++){
            if(inexp[i][j])
            {
                if(inexp[i][j] < 0)
                    e+='!';
                e+=varnames[j];
                if(count-1)
                {
                    e+=K0_D1?" ^ ":" v ";
                    count--;
                }
            }
        }
        e+=')';
        if(e.size()>2){
            exp += e;
            exp += K0_D1?" v ":" ^ ";
        }
    }

    if(exp[exp.size()-2] == 'v' or exp[exp.size()-2] == '^'){
        exp.pop_back();
        exp.pop_back();
    }
}

void boolfunc::MaxRes()
{
    while (Resolution()) {
    }
}

void boolfunc::BlakeAlg(bool KNF0orDNF1)
{
    K0_D1 = KNF0orDNF1;
    if(K0_D1)
        DNF();
    else
        KNF();
    cout<<getExp()<<endl;
    MaxRes();
    cout<<getExp()<<endl;
    MaxAbsorption();
    cout<<getExp()<<endl;
}

string boolfunc::getExp()
{
    return exp;
}

void boolfunc::printINEXP()
{
    int s = inexp.size();
    for(int i = 0; i < s; i++){
        for(int j = 0; j < countvars; j++){
            cout<<inexp[i][j]<<"\t";
        }
        cout<<endl;
    }
    cout<<endl;
}

