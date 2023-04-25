// include file
#include "subspace.h"
#include "deformableParticles2D.h"
#include "DPM_Parallel.h"
#include "random"
#include <cmath>
#include <omp.h>

using namespace std;
extern bool replaceFlag;
extern bool frictionFlag;

// split the packing system into smaller subsystems
void DPM_Parallel::split_into_subspace() {
	int box;
	// create N[0] * N[1] subsystems
	if (subsystem == nullptr)
		subsystem = new frictionlessSubspace[N_systems[0] * N_systems[1]];

	std::vector<double> temp;
	for (int d = 0; d < NDIM; d++)
	{
		double factor = cell(0).pbc.at(d)? 1: 1.2;
		temp.push_back(L.at(d) * factor - BoundaryCoor.at(d));
	}
	// initialize subsystems
	for (int i = 0; i < N_systems[0] * N_systems[1]; i++) {
		subsystem[i].initialize(this, temp, N_systems, i, dt0);
		subsystem[i].cal_cashed_fraction();
		subsystem[i].setFriction(this->coefu, this->coefv);
	}

	// assign vertexes into subsystems
	for (int ci = 0; ci < NCELLS; ci++) {
		for (int vi = 0; vi < cell(ci).getNV(); vi++) {
			cvpair* newpair = new cvpair(ci, vi);
			box = look_for_new_box(newpair);
			migrate_into(box, newpair);
			newpair->boxid = box;
		}
	}
};

// split the packing system into smaller subsystems
void DPM_Parallel_frictionless::split_into_subspace() {
	int box;
	// create N[0] * N[1] subsystems
	//cout << N_systems[0] << ", " << N_systems[1] << endl;
	if (subsystem == nullptr)
		// subsystem = new subspace[N_systems[0] * N_systems[1]];
		subsystem = new frictionlessSubspace[N_systems[0] * N_systems[1]];

	std::vector<double> temp;
	temp.reserve(2);
	for (int d = 0; d < NDIM; d++)
	{
		double factor = cell(0).pbc.at(d)? 1: 1.2;
		temp.push_back(L.at(d) * factor - BoundaryCoor.at(d));
	}
	// initialize subsystems
	for (int i = 0; i < N_systems[0] * N_systems[1]; i++) {
		subsystem[i].initialize(this, temp, N_systems, i, dt0);
		subsystem[i].cal_cashed_fraction();
		subsystem[i].setFriction(this->coefu, this->coefv);
	}

	// assign vertexes into subsystems
	for (int ci = 0; ci < NCELLS; ci++) {
		for (int vi = 0; vi < cell(ci).getNV(); vi++) {
			cvpair* newpair = new cvpair(ci, vi);
			box = look_for_new_box(newpair);
			migrate_into(box, newpair);
			newpair->boxid = box;
		}
	}
};

// cashe list send to subsystems
void DPM_Parallel::cashe_into(int i, vector<cvpair*>& cash_list) {
	subsystem[i].cashe_in(cash_list);
};

// migrate cells into subsystems
void DPM_Parallel::migrate_into(int i, cvpair* const& migration) {
	migration->boxid = i;
	subsystem[i].migrate_in(migration);
};

// figure out which box the cells belong to
int DPM_Parallel::look_for_new_box(cvpair* pair) {
	int box_id = 0;
	int x_id = 0;
	int y_id = 0;

	// figure out x and y index
	x_id = floor(transformPos(pair, 0) / (subsystem[0].L.at(0) / N_systems[0]));
	y_id = floor(transformPos(pair, 1) / (subsystem[0].L.at(1) / N_systems[1]));

	//add periodic boundary just in case
	if (cell(0).pbc.at(0))
		x_id = x_id % N_systems[0];
	// else if (x_id < 0 || x_id >= N_systems[0])
	// 	return -1;
	else if (x_id < 0)
	{
		if (replaceFlag)
			x_id = 0;
		else 
			return -1;
	}
	else if (x_id >= N_systems[0])
	{
		if (replaceFlag)
			x_id = N_systems[0] - 1;
		else 
			return -1;
	}

	if (cell(0).pbc.at(1) || replaceFlag)
		y_id = (y_id + N_systems[1]) % N_systems[1];
	else if (y_id < 0 || y_id >= N_systems[1])
		return -1;

	// convert into box id
	box_id = y_id * N_systems[0] + x_id;

	return box_id;
}

double DPM_Parallel::transformPos(cvpair* pair, int direction) {
	if (cell(pair->ci).pbc.at(direction)) {
		double x = cell(pair->ci).vpos(pair->vi, direction);
		double y = L.at(direction);
		double modPos = x - (int)(x / y) * y;
		return modPos > 0 ? modPos : modPos + y;
	}
	else
	{
		double pos = cell(pair->ci).vpos(pair->vi, direction) - BoundaryCoor.at(direction);
		if (!replaceFlag)
			pos = (!direction && pos < 0)? 1e-10 : pos;
		return pos;
	}
}

// initialization
void DPM_Parallel::initialize_subsystems(int N_x, int N_y) {
	// set how many boxes along each direction
	if (N_systems.size() < 2)
	{
		N_systems.resize(2);
	}

	{
		N_systems.at(0) = N_x;
		N_systems.at(1) = N_y;
	}
	// split
	split_into_subspace();

}


// reset system
void  DPM_Parallel::reset_subsystems() {
	for (int i = 0; i < N_systems.at(0) * N_systems.at(1); i++)
		(subsystem[i]).reset();
}

void  DPM_Parallel::delete_subsystems() {
	delete[] subsystem;
	subsystem = nullptr;
}




void subspace::cal_cashed_fraction() {

	// calculate cashed fraction
	for (int d = 0; d < NDIM; d++) {
		double spacing = L.at(d) / N_systems[d];
		//cashed_fraction.at(d) = pointer_to_system->scale_v(cashed_length) / spacing;
		cashed_fraction.at(d) = 2.5 * pointer_to_system->cell(0).l0 / spacing;
		if (N_systems[d] == 2 && cashed_fraction.at(d) > 0.5) {
			cout << " Too much boxes for too few cells " << endl;
			// has to be exactly 0.5, otherwise there could be problem of untracked bonds
			cashed_fraction.at(d) = 0.5 - 1e-8;
		}
		if (cashed_fraction.at(d) > 1) {
			cout << " Too much boxes for too few cells " << endl;
			cashed_fraction.at(d) = 1 - 1e-8;
		}
		if (this == pointer_to_system->subsystem)
			cout << "fraction at " << d << " = " << cashed_fraction.at(d) << endl;
	}

}


// cashe cells into cashe list
void subspace::cashe_in(vector<cvpair*>& cash_list) {
	cashed_cells.insert(cashed_cells.end(), cash_list.begin(), cash_list.end());
	//for (int i = 0; i < cash_list.size(); i++)
	//	cashed_cells.push_back(cash_list.at(i));
};

// add migrated cells
void subspace::migrate_in(cvpair* const& migration) {
	resident_cells.push_back(migration);
};

// reset cashe system
void subspace::reset_cashe() {
	if (!cashed_cells.empty()) {
		cashed_cells.clear();
	}
}

// reset the whole system
void subspace::reset() {

	reset_cashe();
	if (!resident_cells.empty()) {
		for (auto c : resident_cells)
			delete c;
		resident_cells.clear();
	}
}



// send cashe list 
void subspace::cashe_out(int direction) {

	// find neighbor boxes
	int lower_index = neighbor_box(direction, -1);
	int upper_index = neighbor_box(direction, +1);

	double spacing = L.at(direction) / N_systems[direction];
	// find boundaries of the current box
	double lower_boundary = find_boundary(direction, -1);
	double upper_boundary = find_boundary(direction, +1);

	// reset list
	cash_out_list_lower.clear();
	cash_out_list_upper.clear();

	// check if resident cells are near boundary
	if (!resident_cells.empty()) {
		for (int ci = 0; ci < resident_cells.size(); ci++) {
			if (pointer_to_system->transformPos(resident_cells[ci], direction) < lower_boundary + cashed_fraction.at(direction) * spacing)

				cash_out_list_lower.push_back(resident_cells[ci]);
			// only interact with half of the neighboring subspacese
			//if (resident_cells[ci]->cpos(direction) > upper_boundary - cashed_fraction.at(direction) * spacing)

			//	cash_out_list_upper.push_back(resident_cells[ci]);
		}
	}

	// check if cashed cells are near boundary, but only for y direction
	if (!cashed_cells.empty() && direction == 1) {
		for (int ci = 0; ci < cashed_cells.size(); ci++) {
			if (pointer_to_system->transformPos(cashed_cells[ci], direction) < lower_boundary + cashed_fraction.at(direction) * spacing &&
				pointer_to_system->transformPos(cashed_cells[ci], direction) > lower_boundary)

				cash_out_list_lower.push_back(cashed_cells[ci]);
			if (pointer_to_system->transformPos(cashed_cells[ci], direction) > upper_boundary - cashed_fraction.at(direction) * spacing &&
				pointer_to_system->transformPos(cashed_cells[ci], direction) < upper_boundary)

				cash_out_list_upper.push_back(cashed_cells[ci]);
		}

	}

	// send to other boxes
	if (lower_index >= 0)
		pointer_to_system->cashe_into(lower_index, cash_out_list_lower);
	if (upper_index >= 0)
		pointer_to_system->cashe_into(upper_index, cash_out_list_upper);

};

// find neightor box
int subspace::neighbor_box(int direction, int upper_lower) {
	int current[2], neighbor_box_id;

	// current box id in x and y
	current[0] = box_id % N_systems[0];
	current[1] = floor(box_id / N_systems[0]);

	// neighbor box
	current[direction] += upper_lower;
	// periodic boundary
	if (pointer_to_system->cell(0).pbc.at(direction))
		current[direction] = (current[direction] + N_systems[direction]) % N_systems[direction];
	else if (current[direction] < 0 || current[direction] >= N_systems[direction])
		return -1;
	// neighbor id
	neighbor_box_id = current[0] + current[1] * N_systems[0];

	return neighbor_box_id;
}

double subspace::find_boundary(int direction, int upper_lower) {
	int current[2];
	double boundary;
	double spacing = L.at(direction) / N_systems[direction];

	// current box id in x and y
	current[0] = box_id % N_systems[0];
	current[1] = floor(box_id / N_systems[0]);

	// find boundary
	if (upper_lower == -1)
		boundary = current[direction] * spacing;
	else if (upper_lower == 1)
		boundary = (current[direction] + 1) * spacing;

	return boundary;
}



// migrate cells to other boxes and remove from current box
void subspace::migrate_out() {

	int new_box_index;
	int cell_index;
	cvpair* target_cell;

	// stack indicates near boundary cells that migrate to neighbor boxes
	stack<int> migrate_out_list;
	stack<int> migrate_out_destination;

	/*
	// empty the stack
	while (!migrate_out_list.empty()) {
		migrate_out_list.pop();
	}

	while (!migrate_out_destination.empty()) {
		migrate_out_destination.pop();
	}
	*/

	// find the migration list
	for (int ci = 0; ci < resident_cells.size(); ci++) {
		new_box_index = pointer_to_system->look_for_new_box(resident_cells.at(ci));
		if (new_box_index != box_id) {
			migrate_out_list.push(ci);
			migrate_out_destination.push(new_box_index);
		}
	}


	// migrate to other subsystems
	if (!migrate_out_list.empty()) {
		for (int i = 0; i < migrate_out_list.size(); i++) {
			// migrate backwards, otherwise the list order would be changed
			cell_index = migrate_out_list.top();
			target_cell = resident_cells[cell_index];
			// find which subsystem to go
			new_box_index = migrate_out_destination.top();
			// migrate
			if (new_box_index >= 0)
				pointer_to_system->migrate_into(new_box_index, target_cell);
			else if (replaceFlag)
			{
				cout << "error" << endl;
				cout << target_cell->ci << ", " << target_cell->vi << endl;
				pointer_to_system->printRoutine(0, 10, 0, 0, 0);
				exit(0);
			}
			// pop from list
			migrate_out_list.pop();
			migrate_out_destination.pop();
			// remove from resident list
			resident_cells.erase(resident_cells.begin() + cell_index);
		}
	}
}


// calculate forces 
void subspace::calculateForces_insub() {
	// local variables
	int ci, cj, ck, vi, d, dd, inContact;

	// reset virial stresses to 0
	sigmaXX = 0.0;
	sigmaXY = 0.0;
	sigmaYX = 0.0;
	sigmaYY = 0.0;

	// reset contacts before force calculation
	//resetContacts();
	Ncc = 0;
	Nvv = 0;

	// loop over cells and cell pairs, calculate shape and interaction forces
	if (!resident_cells.empty()) {
		for (ci = 0; ci < resident_cells.size(); ci++) {
			// forces between resident cells
			// loop over pairs, add info to contact matrix
			for (cj = ci + 1; cj < resident_cells.size(); cj++) {
				// if (resident_cells[ci]->ci == resident_cells[cj]->ci)
				// {
				// 	int distance = abs(resident_cells[ci]->vi - resident_cells[cj]->vi);
				// 	if (distance == 1 || distance == (pointer_to_system->cell(resident_cells[ci]->ci).getNV() - 1))
				// 		continue;
				// }
				if (!(pointer_to_system->cell(resident_cells[ci]->ci).inside_hopper) || !(pointer_to_system->cell(resident_cells[cj]->ci).inside_hopper))
					continue;
				if (resident_cells[ci]->boxid != resident_cells[cj]->boxid) {
					cout << "incorrect resident list" << endl;
					//print_information();
					continue;
				}
				// calculate forces, add to number of vertex-vertex contacts
				inContact = vertexForce(resident_cells[ci], resident_cells[cj], sigmaXX, sigmaXY, sigmaYX, sigmaYY);
				if (resident_cells[ci]->ci != resident_cells[cj]->ci && inContact > 0) {
					// add to cell-cell contacts
					if (pointer_to_system->contacts(resident_cells[ci]->ci, resident_cells[cj]->ci) == 0)
					{
						pointer_to_system->addContact(resident_cells[ci]->ci, resident_cells[cj]->ci);
						Ncc++;
					}
					// increment vertex-vertex contacts
					Nvv++;
				}
			}
		}
	}
}

void subspace::calculateForces_betweensub() {
	if (!resident_cells.empty() && !cashed_cells.empty()) {
		for (int ci = 0; ci < resident_cells.size(); ci++) {
			// forces between resident cell and cashed cell
			for (int ck = 0; ck < cashed_cells.size(); ck++) {
				// if (resident_cells[ci]->ci == cashed_cells[ck]->ci)
				// {
				// 	int distance = abs(resident_cells[ci]->vi - cashed_cells[ck]->vi);
				// 	if (distance == 1 || distance == (pointer_to_system->cell(resident_cells[ci]->ci).getNV() - 1))
				// 		continue;
				// }
				if (!(pointer_to_system->cell(resident_cells[ci]->ci).inside_hopper) || !(pointer_to_system->cell(cashed_cells[ck]->ci).inside_hopper))
					continue;
				if (resident_cells[ci]->boxid == cashed_cells[ck]->boxid) {
					cout << "incorrect cashed list" << endl;
					//print_information();
					continue;
				}
				// notice that stress between resident and cashed cells are double counted
				int inContact = vertexForce(resident_cells[ci], cashed_cells[ck], sigmaXX, sigmaXY, sigmaYX, sigmaYY);
				// add to cell-cell contacts
				if (resident_cells[ci]->ci != cashed_cells[ck]->ci && inContact > 0){
					if (pointer_to_system->contacts(resident_cells[ci]->ci, cashed_cells[ck]->ci) == 0)
					{
						pointer_to_system->addContact(resident_cells[ci]->ci, cashed_cells[ck]->ci);
						Ncc++;
					}
					// increment vertex-vertex contacts
					Nvv++;
				}
			}
		}
	}
}
/*
int subspace::vertexForce(cvpair* onTheLeft, cvpair* onTheRight, double& sigmaXX, double& sigmaXY, double& sigmaYX, double& sigmaYY) {
	// return variable
	int inContact = 0;
	deformableParticles2D& leftCell = pointer_to_system->cell(onTheLeft->ci);
	deformableParticles2D& rightCell = pointer_to_system->cell(onTheRight->ci);
	// local variables
	int d, dd;

	// -------------------------
	// 
	// 	   Vertex forces
	//
	// -------------------------

	double forceScale = leftCell.kint;				// force scale
	double energyScale = forceScale;		// energy scale
	double distScale = 0.0;					// distance scale
	double p1 = 1.0 + leftCell.a;					// edge of interaction zone, units of delta
	double ftmp = 0.0;						// temporary force variable
	double uTmp = 0.0;						// temporary energy variable
	double vertexDist = 0.0;				// distance variable
	vector<double> vertexVec(NDIM, 0.0);		// vector to hold vectorial distance quantity
	double contactDistance = 0.0;			// contact distance variable
	double distTmp = 0.0;


	// get distance between vertices i and j
	vertexDist = 0.0;
	for (d = 0; d < NDIM; d++) {
		// get distance to nearest image
		distTmp = leftCell.distance(rightCell, onTheRight->vi, onTheLeft->vi, d);

		// add to vertex distance
		vertexVec.at(d) = distTmp;

		// add to scalar distance
		vertexDist += distTmp * distTmp;
	}

	// get contact distance
	contactDistance = 0.5 * (leftCell.del * leftCell.l0 + rightCell.del * rightCell.l0);

	// get vertex distance
	vertexDist = sqrt(vertexDist);

	// check overlap distances
	if (vertexDist < contactDistance * p1) {
		// increment number of vertex-vertex contacts
		inContact++;

		// define scaled distance (x = distance/contact distance)
		distScale = vertexDist / contactDistance;

		// update force and energy scales
		forceScale = leftCell.kint / contactDistance;
		energyScale = leftCell.kint;

		// IF in zone to use repulsive force (and, if a > 0, bottom of attractive well)
		if (vertexDist < contactDistance) {

			// add to interaction potential
			uTmp = 0.5 * energyScale * pow(1 - distScale, 2) - (energyScale * leftCell.a * leftCell.a) / 6.0;

			leftCell.setUInt(onTheLeft->vi, leftCell.uInt(onTheLeft->vi) + 0.5 * uTmp);
			rightCell.setUInt(onTheRight->vi, rightCell.uInt(onTheRight->vi) + 0.5 * uTmp);

			// add to vectorial forces
			for (d = 0; d < NDIM; d++) {
				// get force value
				ftmp = -forceScale * (1 - distScale) * vertexVec.at(d) / vertexDist;

				// add to force on i
				leftCell.setVForce(onTheLeft->vi, d, leftCell.vforce(onTheLeft->vi, d) + ftmp);

				// subtract off complement from force on j
				rightCell.setVForce(onTheRight->vi, d, rightCell.vforce(onTheRight->vi, d) - ftmp);

				// add to stress tensor
				if (d == 0) {
					sigmaXX -= 2.0 * ftmp * vertexVec.at(0);
					sigmaXY -= 2.0 * ftmp * vertexVec.at(1);
				}
				else {
					sigmaYX -= 2.0 * ftmp * vertexVec.at(0);
					sigmaYY -= 2.0 * ftmp * vertexVec.at(1);
				}
			}
		}

	}

	// return if in contact or not
	return inContact;
}
*/

int subspace::vertexForce(cvpair* onTheLeft, cvpair* onTheRight, double& sigmaXX, double& sigmaXY, double& sigmaYX, double& sigmaYY) {
//int subspace::vertexForce_with_Torque(cvpair* onTheLeft, cvpair* onTheRight, double& sigmaXX, double& sigmaXY, double& sigmaYX, double& sigmaYY) {
	// return variable
	int inContact = 0;
	deformableParticles2D& leftCell = pointer_to_system->cell(onTheLeft->ci);
	deformableParticles2D& rightCell = pointer_to_system->cell(onTheRight->ci);
	// disable this may cause some friction
	// if ((leftCell.vertexEdgeContact[onTheLeft->vi] == onTheRight->ci) && (onTheRight->ci != onTheLeft->ci))
	// 	return 0;
	// local variables
	int d, dd;

	// -------------------------
	// 
	// 	   Vertex forces
	//
	// -------------------------

	double forceScale = leftCell.kint;				// force scale
	double energyScale = forceScale;		// energy scale
	double distScale = 0.0;					// distance scale
	double p1 = 1.0 + leftCell.a;					// edge of interaction zone, units of delta
	double ftmp = 0.0;						// temporary force variable
	double uTmp = 0.0;						// temporary energy variable
	double vertexDist = 0.0;				// distance variable
	vector<double> vertexVec(NDIM, 0.0);		// vector to hold vectorial distance quantity
	double contactDistance = 0.0;			// contact distance variable
	double distTmp = 0.0;


	// get distance between vertices i and j
	vertexDist = 0.0;
	for (d = 0; d < NDIM; d++) {
		// get distance to nearest image
		distTmp = leftCell.distance(rightCell, onTheRight->vi, onTheLeft->vi, d);

		// add to vertex distance
		vertexVec.at(d) = distTmp;

		// add to scalar distance
		vertexDist += distTmp * distTmp;
	}

	// get contact distance
	contactDistance = 0.5 * (leftCell.del * leftCell.l0 + rightCell.del * rightCell.l0) * pointer_to_system->cutoff;

	// get vertex distance
	vertexDist = sqrt(vertexDist);

	// check overlap distances
	if (vertexDist < contactDistance * p1) {
		// increment number of vertex-vertex contacts
		inContact++;

		// define scaled distance (x = distance/contact distance)
		distScale = vertexDist / contactDistance;

		// update force and energy scales
		forceScale = leftCell.kint / contactDistance;
		energyScale = leftCell.kint;

		// IF in zone to use repulsive force (and, if a > 0, bottom of attractive well)
		if (vertexDist < contactDistance) {

			// add to interaction potential
			uTmp = 0.5 * energyScale * pow(1 - distScale, 2) - (energyScale * leftCell.a * leftCell.a) / 6.0;

			leftCell.setUInt(onTheLeft->vi, leftCell.uInt(onTheLeft->vi) + 0.5 * uTmp);
			rightCell.setUInt(onTheRight->vi, rightCell.uInt(onTheRight->vi) + 0.5 * uTmp);

			std::vector<double> force(2);
			std::vector<double> r1(2);
			std::vector<double> r2(2);

			// double normF = forceScale * (1 - distScale) * vertexDist / vertexDist;
			std::vector<double> dv(2);
			std::vector<double> dv_along_norm(2);
			std::vector<double> dv_along_tang(2);
			std::vector<double> friction(2);
			std::vector<double> friction_norm(2);
			std::vector<double> friction_tang(2);
			double projection = 0;
			for (d = 0; d < NDIM; d++) {
				if (frictionFlag){
					dv.at(d) = leftCell.vvel(onTheLeft->vi,d) - rightCell.vvel(onTheRight->vi,d);
					friction.at(d) = - 0.1 * dv.at(d);
					projection += dv.at(d) * vertexVec.at(d) / vertexDist;
				}
				else
				{
					friction.at(d) = 0;
					friction_norm.at(d) = 0;
					friction_tang.at(d) = 0;
				}
			}
			if (frictionFlag){
			// if (0){
				
				for (d = 0; d < NDIM; d++) {
					dv_along_norm.at(d) = projection * vertexVec.at(d) / vertexDist;
					dv_along_tang.at(d) = dv.at(d) - dv_along_norm.at(d);
					friction_norm.at(d) = - coefu * dv_along_norm.at(d);
					friction_tang.at(d) = - coefV * dv_along_tang.at(d);
					// friction_tang.at(d) = 0;
				}
			}
			// cout << "friction force" << endl;
			// cout << friction.at(0) << "," << friction_norm.at(0) + friction_tang.at(0) << endl;
			// add to vectorial forces
			for (d = 0; d < NDIM; d++) {
				// get force value
				ftmp = -forceScale * (1 - distScale) * vertexVec.at(d) / vertexDist;
				force[d] = ftmp;
				r1[d] = leftCell.vrel(onTheLeft->vi, d) + vertexVec.at(d) / 2;
				r2[d] = rightCell.vrel(onTheRight->vi, d) - vertexVec.at(d) / 2;
				// add to force on i
#pragma omp critical
			{
				leftCell.setVForce(onTheLeft->vi, d, leftCell.vforce(onTheLeft->vi, d) + ftmp + friction_norm.at(d) + friction_tang.at(d));

				// subtract off complement from force on j
				rightCell.setVForce(onTheRight->vi, d, rightCell.vforce(onTheRight->vi, d) - ftmp - friction_norm.at(d) - friction_tang.at(d));
			}
				// add to stress tensor
				if (d == 0) {
					sigmaXX -= 2.0 * ftmp * vertexVec.at(0);
					sigmaXY -= 2.0 * ftmp * vertexVec.at(1);
				}
				else {
					sigmaYX -= 2.0 * ftmp * vertexVec.at(0);
					sigmaYY -= 2.0 * ftmp * vertexVec.at(1);
				}
			}
#pragma omp critical
			{
				// torque calculation
				leftCell.torque += r1[0] * force[1] - r1[1] * force[0];
				rightCell.torque += -r2[0] * force[1] + r2[1] * force[0];
			}
		}

	}

	// return if in contact or not
	return inContact;
}

int frictionlessSubspace::vertexEdgeForce(cvpair* onTheLeft, cvpair* onTheRight, double& sigmaXX, double& sigmaXY, double& sigmaYX, double& sigmaYY) {
//int subspace::vertexForce_with_Torque(cvpair* onTheLeft, cvpair* onTheRight, double& sigmaXX, double& sigmaXY, double& sigmaYX, double& sigmaYY) {
	// return variable
	int inContact = 0;
	if (onTheLeft->ci == onTheRight->ci)
	{
		int distance = onTheLeft->vi - onTheRight->vi;
		if (distance == 1 || distance == (-1 * pointer_to_system->cell(onTheRight->ci).getNV() + 1))
			return 0;
	}
	deformableParticles2D& leftCell = pointer_to_system->cell(onTheLeft->ci);
	deformableParticles2D& rightCell = pointer_to_system->cell(onTheRight->ci);
	// local variables
	int d, dd;

	// -------------------------
	// 
	// 	   Vertex forces
	//
	// -------------------------

	double forceScale = leftCell.kint;				// force scale
	double energyScale = forceScale;		// energy scale
	double distScale = 0.0;					// distance scale
	double p1 = 1.0 + leftCell.a;					// edge of interaction zone, units of delta
	double ftmp = 0.0;						// temporary force variable
	double uTmp = 0.0;						// temporary energy variable
	double vertexDist = 0.0;				// distance variable
	vector<double> vertexVec(NDIM, 0.0);		// vector to hold vectorial distance quantity
	double contactDistance = 0.0;			// contact distance variable
	double distTmp = 0.0;


	// get distance between vertices i and j
	vertexDist = vertexEdgeDist(onTheLeft, onTheRight);


	// get contact distance
	contactDistance = 0.5 * (leftCell.del * leftCell.l0 + rightCell.del * rightCell.l0) * pointer_to_system->cutoff;

	if (vertexDist < 0 && onTheLeft->ci != onTheRight->ci)
		return 0;
	// check overlap distances
	// if (vertexDist < contactDistance && vertexDist > -pointer_to_system->insideCutoffFactor * contactDistance) {
	if (vertexDist < contactDistance && vertexDist > -contactDistance) {
		// increment number of vertex-vertex contacts
		inContact++;
		// cvpair onTheLeftV = *onTheLeft;
		// cvpair onTheRightV = *onTheRight;

			// std::vector<cvpair> & onTheRightList = pointer_to_system->collisionMap[onTheLeftV];
			// if (onTheRightList.size() > 0)
			// {
			// 	cout << onTheRightList.size() << endl;
			// 	cout << onTheLeftV.ci << " " << onTheLeftV.vi << " " << onTheLeftV.boxid << endl;
			// 	for (cvpair i : onTheRightList)
			// 	{
			// 		cout << i.ci << " " << i.vi <<  " " << i.boxid <<endl;
			// 	}
			// }
			// onTheRightList.push_back(onTheRightV);
			vector<VECTOR2> v, e;
			v.resize(3);
			VECTOR2 v0, v1, vColission; 
			int nextVi = (onTheRight->vi + 1) % pointer_to_system->cell(onTheRight->ci).getNV();
			for (int d = 0; d < 2; d++)
			{
				vColission[d] = pointer_to_system->cell(onTheLeft->ci).vpos(onTheLeft->vi,d);
				v0[d] = pointer_to_system->cell(onTheRight->ci).vpos(onTheRight->vi,d);
				v1[d] = pointer_to_system->cell(onTheRight->ci).vpos(nextVi,d);
			}
			v[0] = vColission;
			v[1] = v0;
			v[2] = v1;

			e.resize(2);
			e[0] = v[2] - v[1];
			e[1] = v[0] - v[1];
			// deformableParticles2D& leftCell = cell(onTheLeft.ci);
			// deformableParticles2D& rightCell = cell(onTheRight.ci);
			// double eps = 0.5 * (leftCell.del * leftCell.l0 + rightCell.del * rightCell.l0) * cutoff;
			int mode = onTheLeft->ci == onTheRight->ci;
			VECTOR6 forces = -1 * leftCell.kint * pow(contactDistance, -2) * pointer_to_system->gradient(v, e, contactDistance, mode);
			VECTOR6 friction;
			friction.setZero();
			if (frictionFlag){
			// if (0){
				std::vector<double> dv(2);
				std::vector<double> dv_along_norm(2);
				std::vector<double> dv_along_tang(2);
				std::vector<double> frictionT(2);
				std::vector<double> friction_norm(2);
				std::vector<double> friction_tang(2);
				std::vector<double> vAtContact(2);
				VECTOR2 crossed = (-1) * VECTOR2(e[0][1],-e[0][0]);
				double bondL = e[0].norm();
				double distToV1 = e[1].dot(e[0]) / bondL;
				for (d = 0; d < NDIM; d++) {
					vAtContact.at(d) = (rightCell.vvel(onTheRight->vi,d) * (bondL - distToV1) + rightCell.vvel(nextVi,d) * distToV1) / bondL;
					dv.at(d) = leftCell.vvel(onTheLeft->vi,d) - vAtContact.at(d);
					frictionT.at(d) = - 0.1 * dv.at(d);
				}
				double projection = VECTOR2(dv[0],dv[1]).dot(crossed)/crossed.norm();
				for (d = 0; d < NDIM; d++) {
					dv_along_norm.at(d) = projection * crossed[d] / crossed.norm();
					dv_along_tang.at(d) = dv.at(d) - dv_along_norm.at(d);
					friction[d] = - coefu * dv_along_norm.at(d) - coefV * dv_along_tang.at(d);
					friction[d + 2] = - friction[d] * (bondL - distToV1) / bondL;
					friction[d + 4] = - friction[d] * distToV1 / bondL;
				}
				forces += friction;
			}
			// if (onTheLeft->ci == 0 || onTheRight->ci == 0)
			// 	if (forces.maxCoeff() > 1)
			// 	{
			// 		cout << "large force" << endl;
			// 		cout << onTheLeft->ci << " " << onTheRight->ci << endl;
			// 			for (int i = 0; i < 6; i++)
			// 				cout << forces[i];
			// 	}
			// VECTOR6 absF = forces.cwiseAbs();
			// double maxF = absF.maxCoeff();
			// if (maxF > 1e3)
			// {
			// 	// forces /= maxF;
			// 	cout << "force error" << endl;
			// 	cout << onTheLeft->ci << ", " << onTheLeft->vi << ", " << onTheRight->ci << ", " << onTheRight->vi << endl;
			// 	#pragma omp critical
			// 	{
			// 		pointer_to_system->printRoutine(0, 10, 0, 0, 0);
			// 	}
			// 	exit(0);
			// }
		#pragma omp critical
		{
			leftCell.vertexEdgeContact[onTheLeft->vi] = onTheRight->ci;
			for (int d = 0; d < 2; d++)
			{
				leftCell.setVForce(onTheLeft->vi, d, leftCell.vforce(onTheLeft->vi, d) + forces[d]);
				rightCell.setVForce(onTheRight->vi, d, rightCell.vforce(onTheRight->vi, d) + forces[d + 2]);
				rightCell.setVForce(nextVi, d, rightCell.vforce(nextVi, d) + forces[d + 4]);
			}
		}
	}

	// return if in contact or not
	return inContact;
}

double frictionlessSubspace::vertexEdgeDist(const cvpair* onTheLeft, const cvpair* onTheRight)
{
	VECTOR2 v0, v1, v; 
	for (int d = 0; d < 2; d++)
	{
		v[d] = pointer_to_system->cell(onTheLeft->ci).vpos(onTheLeft->vi,d);
		v0[d] = pointer_to_system->cell(onTheRight->ci).vpos(onTheRight->vi,d);
		int nextVi = (onTheRight->vi + 1) % pointer_to_system->cell(onTheRight->ci).getNV();
		v1[d] = pointer_to_system->cell(onTheRight->ci).vpos(nextVi,d);
	}
	
  const VECTOR2 e = v1 - v0;
  const VECTOR2 e1 = v - v0;
  const VECTOR2 n = VECTOR2(e[1],-e[0]);
 
  const double check = e1.dot(e);

  // if the point projects to inside the segment
  if (check > 0 && check < e.dot(e))
  {
    const VECTOR2 nHat = n / n.norm();
    const double normalDistance = (nHat.dot(v - v0));
    return normalDistance;
  }
  else
	return 1e8;

  // get the distance to each vertex
  const VECTOR2 vertexDistances((v - v0).norm(), 
                                (v - v1).norm());

  // get the smallest of both the edge and vertex distances
  const double vertexMin = vertexDistances.minCoeff();

  // return the smallest of those
  return vertexMin;
}


// calculate forces 
void frictionlessSubspace::calculateEdgeForces_insub() {
	// local variables
	int ci, cj, ck, vi, d, dd, inContact;

	// reset virial stresses to 0
	sigmaXX = 0.0;
	sigmaXY = 0.0;
	sigmaYX = 0.0;
	sigmaYY = 0.0;

	// reset contacts before force calculation
	//resetContacts();
	Ncc = 0;
	Nvv = 0;

	// loop over cells and cell pairs, calculate shape and interaction forces
	if (!resident_cells.empty()) {
		for (ci = 0; ci < resident_cells.size(); ci++) {
			// forces between resident cells
			// loop over pairs, add info to contact matrix
			// int edgeContactForCi = 0;
			for (cj = 0; cj < resident_cells.size(); cj++) {
				if (ci == cj) continue;
				// if (resident_cells[ci]->ci == resident_cells[cj]->ci)
				// {
				// 	int distance = resident_cells[ci]->vi - resident_cells[cj]->vi;
				// 	if (distance == 1 || distance == (-1 * pointer_to_system->cell(resident_cells[ci]->ci).getNV() + 1))
				// 		continue;
				// }
				if (!(pointer_to_system->cell(resident_cells[ci]->ci).inside_hopper) || !(pointer_to_system->cell(resident_cells[cj]->ci).inside_hopper))
					continue;
				if (settleDown && (!(pointer_to_system->cell(resident_cells[ci]->ci).inside_hopper==2) && !(pointer_to_system->cell(resident_cells[cj]->ci).inside_hopper==2)))
					continue;
				if (resident_cells[ci]->boxid != resident_cells[cj]->boxid) {
					cout << "incorrect resident list" << endl;
					//print_information();
					continue;
				}
				// calculate forces, add to number of vertex-vertex contacts
				inContact = vertexEdgeForce(resident_cells[ci], resident_cells[cj], sigmaXX, sigmaXY, sigmaYX, sigmaYY);
				// edgeContactForCi += inContact;
				if (resident_cells[ci]->ci != resident_cells[cj]->ci && inContact > 0) {
					// add to cell-cell contacts
					if (pointer_to_system->contacts(resident_cells[ci]->ci, resident_cells[cj]->ci) == 0)
					{
						pointer_to_system->addContact(resident_cells[ci]->ci, resident_cells[cj]->ci);
						Ncc++;
					}
					// increment vertex-vertex contacts
					Nvv++;
				}
			}
		}
	}
}

void frictionlessSubspace::calculateEdgeForces_betweensub() {
	if (!resident_cells.empty() && !cashed_cells.empty()) {
		for (int ci = 0; ci < resident_cells.size(); ci++) {
			// forces between resident cell and cashed cell
			for (int ck = 0; ck < cashed_cells.size(); ck++) {
				// if (resident_cells[ci]->ci == cashed_cells[ck]->ci)
				// {
				// 	int distance = resident_cells[ci]->vi - cashed_cells[ck]->vi;
				// 	if (distance == 1 || distance == (-1 * pointer_to_system->cell(resident_cells[ci]->ci).getNV() + 1))
				// 		continue;
				// }
				if (!(pointer_to_system->cell(resident_cells[ci]->ci).inside_hopper) || !(pointer_to_system->cell(cashed_cells[ck]->ci).inside_hopper))
					continue;
				if (settleDown && (!(pointer_to_system->cell(resident_cells[ci]->ci).inside_hopper==2) && !(pointer_to_system->cell(cashed_cells[ck]->ci).inside_hopper==2)))
					continue;
				if (resident_cells[ci]->boxid == cashed_cells[ck]->boxid) {
					cout << "incorrect cashed list" << endl;
					//print_information();
					continue;
				}
				// notice that stress between resident and cashed cells are double counted
				int inContact = vertexEdgeForce(resident_cells[ci], cashed_cells[ck], sigmaXX, sigmaXY, sigmaYX, sigmaYY);
				inContact += vertexEdgeForce(cashed_cells[ck], resident_cells[ci], sigmaXX, sigmaXY, sigmaYX, sigmaYY);
				// add to cell-cell contacts
				if (resident_cells[ci]->ci != cashed_cells[ck]->ci && inContact > 0){
					if (pointer_to_system->contacts(resident_cells[ci]->ci, cashed_cells[ck]->ci) == 0)
					{
						pointer_to_system->addContact(resident_cells[ci]->ci, cashed_cells[ck]->ci);
						Ncc++;
					}
					// increment vertex-vertex contacts
					Nvv++;
				}
			}
		}
	}
}
