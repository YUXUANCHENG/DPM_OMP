#ifndef TAOSOLVER_H
#define TAOSOLVER_H

#include "cellPacking2D.h"
#include "NVE.h"

class TaoSolver {
public:

	cellPacking2D* cellpointer = nullptr;
	DPMNVEsimulator* simulator = nullptr;

	TaoSolver(cellPacking2D* cell, DPMNVEsimulator* sim) {
		cellpointer = cell;
		simulator = sim;
	}

	double* NVE_tao(double T, double v0, double Dr, double vtau, double t_scale, int frames) {
		cellpointer->print_frequency = ceil(T / (cellpointer->dt0 * t_scale * frames));

		//sp_NVE(T, v0, Dr, vtau, t_scale, frames);
		double scaled_v = cellpointer->scale_v(v0);
		double v_x = 0, v_y = 0, cur_speed = 0;
		for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
			v_x = cellpointer->cell(ci).cal_mean_v(0);
			v_y = cellpointer->cell(ci).cal_mean_v(1);
			cur_speed += sqrt(v_x * v_x + v_y * v_y);
		}
		cur_speed /= cellpointer->NCELLS;
		double factor = sqrt(2) * scaled_v / cur_speed;
		cellpointer->rescaleVelocities(factor * factor * cellpointer->totalKineticEnergy());

		int Ntotal = floor(T / (cellpointer->dt0 * t_scale));

		std::vector< std::vector<double>> x_com(frames, std::vector<double>(cellpointer->NCELLS, 0));
		std::vector< std::vector<double>> y_com(frames, std::vector<double>(cellpointer->NCELLS, 0));
		std::vector<double> speed(frames, 0);
		double* results = new double[2];
		bool REACHED = false;
		bool READY = false;
		std::vector<double> tao;
		std::vector<double> Temp;
		int count = 0;

		while (true)
		{
			int rowIndex = 0;
			int colIndex = 0;
			int taos[2];
			double q = PI / sqrt(cellpointer->L[0] * cellpointer->L[1] * cellpointer->phi / (PI * cellpointer->NCELLS));
			for (int r = 0; r < 2; r++)
			{
				rowIndex = 0;
				for (int count = 0; count < Ntotal; count++)
				{
					simulator->NVERoutine();

					if (count % cellpointer->print_frequency == 0 && rowIndex < frames)
					{
						for (int ci = 0; ci < cellpointer->NCELLS; ci++) {
							v_x = cellpointer->cell(ci).cal_mean_v(0);
							v_y = cellpointer->cell(ci).cal_mean_v(1);
							cur_speed += sqrt(v_x * v_x + v_y * v_y);
						}
						cur_speed /= cellpointer->NCELLS;
						speed[rowIndex] = cur_speed;

						for (int i = 0; i < cellpointer->NCELLS; i++)
						{
							x_com[rowIndex][i] = cellpointer->cell(i).cpos(0);
							y_com[rowIndex][i] = cellpointer->cell(i).cpos(1);
						}
						rowIndex++;
					}
				}
				taos[r] = cellpointer->calTao(q, rowIndex, x_com, y_com);
			}
			double meanSpeed = 0;
			for (int m = 0; m < rowIndex; m++)
				meanSpeed += speed[m];
			meanSpeed /= rowIndex;
			cout << omp_get_thread_num() << " : current tao: " << taos[0] * cellpointer->dt0 * cellpointer->print_frequency << "," << taos[1] * cellpointer->dt0 * cellpointer->print_frequency << endl;
			if (taos[0] > 0 && taos[1] > 0)
			{
				//double testEq = double(abs(taos[1] - taos[0])) / taos[0];
				double testEq = abs((log(taos[1]) - log(taos[0])) / log(taos[0] * cellpointer->dt0 * cellpointer->print_frequency * meanSpeed));
				if ((taos[0] + taos[1]) / 2 > exp(log(Ntotal) * 2.0 / 3.0))
				{
					Ntotal *= 5;
					cellpointer->print_frequency *= 5;
					cout << omp_get_thread_num() << " : T is a bit small" << endl;
				}
				else if ((taos[0] + taos[1]) / 2 < exp(log(Ntotal) / 3.0))
				{
					Ntotal /= 5;
					cellpointer->print_frequency /= 5;
					cout << omp_get_thread_num() << " : T is a bit large" << endl;
				}
				else if ((testEq < 0.05 && !REACHED) || REACHED)
					//else if (testEq < 0.5)
				{
					//if (count > 1){
					tao.push_back(cellpointer->dt0 * cellpointer->print_frequency * (taos[1] + taos[0]) / 2);
					Temp.push_back(meanSpeed);
					//}
					if (!REACHED)
					{
						cout << omp_get_thread_num() << " : Equilibrated" << endl;
						REACHED = true;
					}
					else
						cout << omp_get_thread_num() << " : Extended runs" << endl;
					count++;
					if (count > 3) READY = true;
					//if (count > 5) READY = true;
				}
				if (READY)
				{
					double avg = 0;
					for (double item : tao)
						avg += item;
					avg /= tao.size();
					results[0] = avg;
					avg = 0;
					for (double item : Temp)
						avg += item;
					avg /= Temp.size();
					results[1] = avg;
					return results;
				}

			}
			//else if(taos[0] < 0 && taos[1] < 0){
			else {
				Ntotal *= 5;
				cellpointer->print_frequency *= 5;
				cout << omp_get_thread_num() << " : T is too small" << endl;
			}
		}
	}


};

#endif