#include "cell_jamming.h"
#include "CLI.h"
#include "CLI_hopper.h"
#include "DPM_Parallel.h"
#include "bumpy_Parallel.h"
#include "bumpyEllipse_Parallel.h"
#include "bumpyDimer_Parallel.h"



// use std name space
using namespace std;
extern bool pForceFlag;
bool pForceFlag = 1;

int main(int argc, char const *argv[]){

# pragma omp parallel
{
	# pragma omp master
	cout << "# of threads = " << omp_get_num_threads() << endl;
}

	//jamming main_function;
	//main_function.sp_NVE_arg(argv);
	//main_function.sp_NVE_tao(argv);

	//Bumpy_CLI<> cli;
	//Bumpy_CLI<Bumpy_Parallel> cli;
	//BumpyDimer_CLI<> cli;
	//BumpyDimer_CLI<BumpyDimer_Parallel> cli;
	//DPM_CLI<> cli;
	DPM_CLI<DPM_Parallel> cli;
	//BumpyEllipse_CLI<> cli;
	//BumpyEllipse_CLI<BumpyEllipse_Parallel> cli;
	
	//DPM_Hopper_CLI<> cli;
	//DPM_Hopper_CLI<DPM_Parallel> cli;
	//Bumpy_Hopper_CLI<> cli;
	//Bumpy_Hopper_CLI<BumpyDimer_Parallel> cli;

	//cli.findJamming(argv);
	//cli.NVE(argv);
	cli.NVEvsDPhi(argv);
	//cli.NVEvsDPhi(argv);

	//cli.hopperFlow(argv);
	//cli.deformation(argv);
	
	//cli.calTao(argv);
	//cli.ArrheniusAngell(argv);

	system("pause");
	return 0;

}


