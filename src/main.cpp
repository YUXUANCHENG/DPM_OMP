#include "cell_jamming.h"
#include "CLI.h"
#include "CLI_hopper.h"
#include "DPM_Parallel.h"
#include "bumpy_Parallel.h"
#include "bumpyEllipse_Parallel.h"
#include "bumpyDimer_Parallel.h"



// use std name space
using namespace std;
extern bool wallAttracFlag;
bool wallAttracFlag = 0;
extern bool constPressureFlag;
bool constPressureFlag = 0;

extern bool replaceFlag;
bool replaceFlag = 1;
extern bool frictionFlag;
bool frictionFlag = 1;

extern bool variableExtFflag;
bool variableExtFflag = 0;
extern bool horrizontalPistonFlag;
bool horrizontalPistonFlag = 0;
extern bool frictionalWallFlag;
bool frictionalWallFlag = 0;

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
	//DPM_CLI<DPM_Parallel> cli;
	//BumpyEllipse_CLI<> cli;
	//BumpyEllipse_CLI<BumpyEllipse_Parallel> cli;
	
	// DPM_Hopper_CLI<> cli;
	DPM_Hopper_CLI<DPM_Parallel> cli;
	// Bumpy_Hopper_CLI<> cli;
	// Bumpy_Hopper_CLI<Bumpy_Parallel> cli;

	//cli.findJamming(argv);
	//cli.NVE(argv);
	//cli.NVEvsDPhi(argv);

	cli.hopperFlow(argv);
	// cli.measureFriction(argv);
	// cli.deformation(argv);
	//cli.findParameter(argv);
	
	//cli.calTao(argv);
	//cli.ArrheniusAngell(argv);

	//optFun();

	//system("pause");
	return 0;

}


