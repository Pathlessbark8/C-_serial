#include<fstream>
#include<random>

using namespace std;

int main(){
    srand(time(nullptr));
    int MAX=3;
    int MIN=0;
    fstream fout;
    fout.open("data.dat",ios::out);
    for(int i=0;i<10000;i++)
    {
        fout<<MIN+(double)(rand()) / ((double)(RAND_MAX/(MAX - MIN)))<<" "<<MIN+(double)(rand()) / ((double)(RAND_MAX/(MAX - MIN)))<<" "<<MIN+(double)(rand()) / ((double)(RAND_MAX/(MAX - MIN)))<<" "<<MIN+(double)(rand()) / ((double)(RAND_MAX/(MAX - MIN)))<<" ";
        fout<<MIN+(double)(rand()) / ((double)(RAND_MAX/(MAX - MIN)))<<" ";
        fout<<MIN+(double)(rand()) / ((double)(RAND_MAX/(MAX - MIN)))<<" ";
        fout<<MIN+(double)(rand()) / ((double)(RAND_MAX/(MAX - MIN)))<<" ";
        fout<<MIN+(double)(rand()) / ((double)(RAND_MAX/(MAX - MIN)))<<" ";
        fout<<MIN+(double)(rand()) / ((double)(RAND_MAX/(MAX - MIN)))<<" ";
        fout<<MIN+(double)(rand()) / ((double)(RAND_MAX/(MAX - MIN)))<<"\n";
    }
    fout.close();
}