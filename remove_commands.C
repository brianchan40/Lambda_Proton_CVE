#include <fstream>
#include <iostream>
#include <TFile>

ofstream result_file;

void remove_commands(int option){
    
    system("rm running.txt");
    system("touch running.txt");
    if(option == 0) system("condor_q brian40 | grep R >> running.txt");
    else if(option == 1) system("condor_q brian40 | grep 74.0 >> running.txt");    

    result_file.open("remove_commands.sh");
    std::fstream eff_lam_file("./running.txt", std::ios_base::in);
    
    string a;
    int count = 0;
    while(eff_lam_file >> a){
        if(option == 0) {if((count - 8)%12 == 0) {result_file << "condor_rm " << a << "\n";}}
	else if(option == 1) {if((count)%12 == 0) {result_file << "condor_rm " << a << "\n";}}
//        cout << a << endl;
        count++;
    }
    
    result_file.close();
    
    system("chmod 777 ./remove_commands.sh");
    system("./remove_commands.sh");
}
