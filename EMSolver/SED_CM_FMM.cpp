#include "SED_CM_FMM.h"
#include "ConfigLoader.h"
#include "Mathdef.h"
#include "Quadrature.h"
#include "tools.h"
#include "CommonEdge.h"
#include "CommonTriangle.h"
#include "Triangle.h"
#include "Tetrahedron.h"
#include "iml.h"
#include "Miscellaneous.h"


using namespace component;
using namespace mom;
using namespace math;
using std::setw;
using Eigen::SparseMatrixBase;


SED_CM_FMM::SED_CM_FMM()
{

}

SED_CM_FMM::~SED_CM_FMM()
{
}

void SED_CM_FMM::init(component::ConfigLoader * ploader)
{
	EM::init(ploader);
	ex_num = 0;

	if (rtype_ == policy::ResultType::Rad)
		readExcEdges(dir_ + '/' + folder_name_ + ".rad");

	dir_ = tool::creatFolder(dir_, folder_name_ + "_SED_CM_FMM");

	//FEKOcur_path = ploader->getFEKOcurPath();
	//auto pos_suffix = FEKOcur_path.find_last_of('.');
	//if (pos_suffix == Qstring::npos)
		//throw ConfigInvalid("FEKO current path parameter cannot match");

	auto log_name = dir_ + (rtype_ == policy::ResultType::Rad ? "/rad_" : "/sca_") + "runtime.log";
	logger_.open(log_name);

	mesh_ptr_ = std::make_shared<Mesh>();
	ce_ptr_ = std::make_shared<CommonEdge>();
	ct_ptr_ = std::make_shared<CommonTriangle>();

	isContinuous = ploader->getisContinuous();
	Sigma = ploader->getSigma();
	Dx = ploader->getDx();
	Dy = ploader->getDy();
	Array_x = ploader->getArray_x();
	Array_y = ploader->getArray_y();
	Array_num = Array_x * Array_y;
	Phase_0 = ploader->getPhase_0();
	Phase_x = ploader->getPhase_x();
	Phase_y = ploader->getPhase_y();
	scale_x = ploader->getScale_alongx();
	/////////////////////////////////
	k_ = PI2 * incidence_.freq / cc;
	omiga = PI2 * incidence_.freq;
	/////////////////////////////////

	/////////////////////////////////
	max_iter_num_ = ploader->getMaxIterationNum();
	iter_threshold = ploader->getIterationThreshold();
	svd_threshold = ploader->getACAThreshold();
	preconditioning_ = ploader->getPreconditioningSwitch();
	/////////////////////////////////

	////////////////////////////////
	multiInc = ploader->getMultiInc();
	if (multiInc.PW_theta_to < multiInc.PW_theta_from)
		multiInc.PW_num = 0;
	else
	{
		int theta_num = ((multiInc.PW_theta_to - multiInc.PW_theta_from) / multiInc.PW_theta_delta) + 1;
		int phi_num = ((multiInc.PW_phi_to - multiInc.PW_phi_from) / multiInc.PW_phi_delta) + 1;
		multiInc.PW_num = (multiInc.polarization == "Bidirect") ? 2 * theta_num*phi_num : theta_num * phi_num;
	}
	////////////////////////////////

	mesh_ptr_->loadVSIEMesh(ploader->getMeshPath());
	mesh_ptr_->getBoundary(min_boundary, max_boundary);
	cen_box = mesh_ptr_->getCenterBox();
	//feko_mesh_ptr_->loadMeshFile(FEKOcur_path.substr(0, FEKOcur_path.find_last_of('.')) + ".nas");
	/////////////////////////////////////
	ce_ptr_->buildSEDCommonEdge_com(mesh_ptr_);
	ct_ptr_->buildSEDCommonTriangle_com(mesh_ptr_, isContinuous);

	threshold_edm = 1.0f*cc / incidence_.freq;
	bou_ce1 = 0;
	bou_ce2 = 0;
	if (isContinuous == 1)
	{
		ce_ptr_->buildSEDCommonEdge_boundary(mesh_ptr_, Dx, Dy);
		ce_ptr_->getBoundary_ce(bou_ce1, bou_ce2);
		ct_ptr_->buildSEDCommonTriangle_boundary(mesh_ptr_, Dx, Dy);
	}

	ce_ptr_->combineSEDCommonEdge(mesh_ptr_);
	ct_ptr_->combineSEDCommonTriangle(mesh_ptr_);

	com_ce = ce_ptr_->getCommon_ce();
	com_ct = ct_ptr_->getCommon_ct();
	///////////////////////////////////////////

	unknowns_t = ct_ptr_->getCommonTriangleNum();
	unknowns_e = ce_ptr_->getCommonEdgeNum();
	unknowns_ = unknowns_t + unknowns_e;
	tet_num = mesh_ptr_->getTetrahedronNum();
	///////////////////////////////////////////////////////

	auto kd = k_ * sqrt(Dx*Dx + Dy * Dy);
	sam_theta_ = static_cast<int>(kd + 6 * cbrtf(kd));
	sam_phi_ = 2 * sam_theta_;

	////////////////////////////////Paralle/////////////////////////////////////
	pool_.init(ploader->getThreadNumber());
	auto task_num = pool_.getThreadNum() * ploader->getTaskFactor();

	Task_construct();
	Task_assign(task_num);
	
	////////////////////////////////////////////////////////////////////////////
	logger_ << ploader->getReport();
	reportInfo(logger_);
	TIME_LOG("init");
}

void SED_CM_FMM::prepareFSPGF(value_t _Dx, value_t _Dy, int _Array_x, int _Array_y, int _t_sum)
{
	auto lamda = cc / incidence_.freq;
	Dx = _Dx * lamda;
	Dy = _Dy * lamda;
	Array_x = _Array_x;
	Array_y = _Array_y;
	t_sum = _t_sum;
}

void SED_CM_FMM::solve()
{
	SEGMENT("Solve");
	//Z.zeros(unknowns_, unknowns_);
	preCalculateKArray();

	Z_sed.zeros(unknowns_, unknowns_);
	Z_sed_eigen.resize(unknowns_, unknowns_);

	I_sed.zeros(unknowns_, multiInc.PW_num + ex_num);

	//I_sed_eigen_rad.resize(unknowns_, multiInc.PW_num + 1);

	V.zeros(unknowns_);
	V_sed.zeros(unknowns_, multiInc.PW_num);
	I_feed.resize(Array_x*Array_y);

	
	//readFEKOcurrent();
	//multiInc.PW_num = 0;
	if (rtype_ == policy::ResultType::Rad)
	{
		I_sed_eigen.resize(unknowns_, multiInc.PW_num + ex_num);
		V_sed_eigen.setZero(unknowns_, multiInc.PW_num + ex_num);
	}
	else
	{
		I_sed_eigen.resize(unknowns_, multiInc.PW_num);
		V_sed_eigen.resize(unknowns_, multiInc.PW_num);
	}

	Z2K1.zeros(tet_num, tet_num);
	Z2K2.set_size(tet_num, tet_num);
	Z2K3.set_size(tet_num, tet_num);
	Z2K4.zeros(tet_num, tet_num);
	Z5K.zeros(tet_num, tet_num);

	RIVP.zeros(tet_num);
	RRIV.zeros(tet_num);
	IV.zeros(tet_num);

	RIV.resize(tet_num);
	IVP.resize(tet_num);

	/////////////SED/////////////////////
	LOG(fillSEDZ_MP(), "Filling 3x3 array impedance matrix");//yigai
	TIME_LOG("fillSEDZ");

	if (multiInc.PW_num != 0)
		LOG(fillSEDV(), "Filling 3x3 r.h.s");//yigai

	if (rtype_ == policy::ResultType::Rad&&ex_num != 0)
	{
		LOG(fillSEDV_rad(), "Filling radiation r.h.s");
	}
	TIME_LOG("fillSEDV");

	//LOG(getJCM(), "Construct current by CM");
	//TIME_LOG("Construct the current by CM");
	LOG(solveIsed(), "Construct current of 3x3 array");
	TIME_LOG("solveIsed");

	//////////////REDUCTION/////////////////
	LOG(fillSVDI(), "Filling the current post SVD");
	TIME_LOG("fillSVDI");

	LOG(buildblk(), "Filling the BlockArray");
	TIME_LOG("buildblk");

	LOG(preCalculateNearZ(), "preCalculating near interaction");
	TIME_LOG("preCalculateNearZ");

	LOG(CalculateNearZ_v2(), "Calculating near interaction");
	TIME_LOG("CalculateNearZ");

	if (preconditioning_)
	{
		LOG(calculatePreconditioner(), "Calculating preconditioner");
		TIME_LOG("Preconditioner");
	}

	LOG(preCalculateRadAndRecv_v2(), "Calculating rad and recv");
	TIME_LOG("CalculateRadAndRecv");

	LOG(preCalculateTrans_v2(), "Calculating transfer");
	TIME_LOG("CalculateTransfer");

	//LOG(fillZbyFMM(), "Calculating ZbyFMM");
	//TIME_LOG("fillZbyFMM");
	//LOG(fillREDZ_MP(), "Filling the reduced matrix");
	//LOG(fillREDZ(), "Filling the reduced matrix");//yigai
	//TIME_LOG("fillREDZ");

	if (rtype_ == policy::ResultType::Rad)
	{
		LOG(fillREDV_rad(), "Filling the reduced r.h.s");
	}
	else
	{
		LOG(fillREDV_MP(), "Filling the reduced r.h.s");//yigai
	}
	TIME_LOG("fillREDV");

	Qcout << setw(30) << "Iterative solving:" << std::endl;
	int iter_num = iterativeSolve_GMRES();
	if (iter_num == max_iter_num_)
		logger_ << LEVEL1 "WARNING: Iteration cannot converge" << std::endl;
	TIME_LOG("iterativeSolve");

	int useless = 0;
	//if (!iml::BICGSTAB(Z_red, V_red, I_red, iter_threshold, iter_num, useless))
	//throw std::runtime_error("fail solving matrix equation");
	//if (!arma::solve(I_red, Z_red, V_red, arma::solve_opts::fast))
	//throw std::runtime_error("fail solving matrix equation");
	//Qcout << "Begin the solving the reduced matrix..." << std::endl;//////////////////////
	//I_red_eigen = Z_red_eigen.partialPivLu().solve(V_red_eigen);////////////////////////////
	//Eigen::GMRES<CMatrix> bicg;
	//bicg.setTolerance(0.001);
	//bicg.compute(Z_red_eigen);
	//I_red_eigen = bicg.solve(V_red_eigen);
	//if(!iml::BICGSTAB(Z_red_eigen, V_red_eigen, I_red_eigen, iter_threshold,500, iter_num))
	//throw std::runtime_error("fail solving matrix equation");
	Qcout << "Finish the solving it!" << std::endl;
	TIME_LOG("solve");

	LOG(constructI(), "Construct the finally current...");
	Qcout << "Finish the constructing!" << std::endl;
	TIME_LOG("cosntructI");

	////////////////////////////////////////////////////////////
	//I.zeros(unknowns_);


	/*if (rtype_ == policy::ResultType::Rad)
	{
	LOG(radiateV(), "Filling rad voltage vector");
	}
	else
	{
	LOG(fillV(), "Filling sca voltage vector");
	}
	TIME_LOG("fillV");*/

	//if (!arma::solve(I, Z, V, arma::solve_opts::fast))
	//if (!iml::BICGSTAB(Z, V, I, iter_threshold, iter_num, useless))
	//throw std::runtime_error("fail solving matrix equation");
	Qcout << "Saving the data..." << std::endl;

	if (writeZIVData())
		Qcout << "The data is already saved!" << std::endl;
	Qcout << "success" << std::endl;
	TIME_LOG("Save data");


#ifdef _DEBUG
	//if (!writeZIVData())
	//  Qcerr << "fail writing ZIV data" << std::endl;
#endif
}

void SED_CM_FMM::output()
{
	SEGMENT("Result");
	Result result;
	if (rtype_ == policy::ResultType::Rad)
	{
		LOG(result.getRadiationField(this, &SED_CM_FMM::getEFiled), "Calculating Radiation");
		TIME_LOG("getRadiation");
		LOG(getCurrentOnFeed(), "Calculating the current on the feed point");
	}
	else
	{
		LOG(result.getBistaticRCS(this, &SED_CM_FMM::getBiRCS), "Calculating BistaticRCS");
		TIME_LOG("getBistaticRCS");
	}
	LOG(result.getNearField(this, &SED_CM_FMM::getNearEField), "Calculating Near field");

	LOG(result.getCurrentDistribution(this, &SED_CM_FMM::calculateSurfaceCurrent), "Calculating surface current");
}

void SED_CM_FMM::clear()
{
	pool_.clear();
	Z.reset();
	I.reset();
	V.reset();
	exc_edge_.clear();

	mesh_ptr_->clear();
	ct_ptr_->clear();
}

void SED_CM_FMM::reportInfo(Qostream& strm) const
{
	mesh_ptr_->reportInfo(strm);
	ct_ptr_->reportInfo(strm);
}

void SED_CM_FMM::Task_construct()
{
	for (int f = 0; f < 9; f++)
	{
		for (int s = 0; s < 9; s++)
		{
			std::pair<int, int> task_temp(f, s);
			SED_task1.push_back(task_temp);
		}
	}

	for (int f = 1 - Array_x; f < Array_x; f++)
	{
		for (int s = 1 - Array_y; s < Array_y; s++)
		{
			std::pair<int, int> task_temp(f, s);
			SED_task2.push_back(task_temp);
		}
	}

	int array_num = Array_x * Array_y;
	for (int v = 0; v < array_num; v++)
	{
		REDV_task.push_back(v);
	}

	for(int idx=0;idx<Array_x;idx++)
		for (int idy = 0; idy < Array_y; idy++)
		{
			std::pair<int, int> task_temp(idx, idy);
			rad_rec_task.push_back(task_temp);
		}

	for (int f = 0; f < array_num; f++)
	{
		int idxm = f / Array_y;//test
		int idym = f % Array_y;//test
		for (int s = 0; s < array_num; s++)
		{
			int idxn = s / Array_y;//test
			int idyn = s % Array_y;//test
			VectorR3 dis(Dx*(idxm - idxn), Dy*(idym - idyn), 0.0);
			if (dis.Norm() < 9.0*lamda)
			{
				std::pair<int, int> task_temp(f, s);
				trans_task.push_back(task_temp);
			}
			else
			{
				std::pair<int, int> task_temp(f, s);
				dp_task.push_back(task_temp);
			}
			
		}
	}
		
}

void SED_CM_FMM::Task_assign(size_t task_num)
{
	auto sed1_pair_num = SED_task1.size();
	auto sed2_pair_num = SED_task2.size();
	auto redv_num = Array_x * Array_y;
	auto rad_rec_num = rad_rec_task.size();
	auto trans_num = trans_task.size();
	auto dp_num = dp_task.size();

	int average1 = sed1_pair_num / task_num;
	int remaining1 = sed1_pair_num % task_num;
	int cur_num1 = 0;

	int average2 = sed2_pair_num / task_num;
	int remaining2 = sed2_pair_num % task_num;
	int cur_num2 = 0;

	int average3 = redv_num / task_num;
	int remaining3 = redv_num % task_num;
	int cur_num3 = 0;

	int average4 = rad_rec_num / task_num;
	int remaining4 = rad_rec_num % task_num;
	int cur_num4 = 0;

	int average5 = trans_num / task_num;
	int remaining5 = trans_num % task_num;
	int cur_num5 = 0;

	int average6 = dp_num / task_num;
	int remaining6 = dp_num % task_num;
	int cur_num6 = 0;

	for (int i = 0; i <= task_num; i++)
	{
		Task_assign1.push_back(cur_num1);
		cur_num1 += average1;
		if (remaining1)
		{
			++cur_num1;
			--remaining1;
		}

		Task_assign2.push_back(cur_num2);
		cur_num2 += average2;
		if (remaining2)
		{
			++cur_num2;
			--remaining2;
		}

		Task_assign3.push_back(cur_num3);
		cur_num3 += average3;
		if (remaining3)
		{
			++cur_num3;
			--remaining3;
		}

		rad_rec_assign.push_back(cur_num4);
		cur_num4 += average4;
		if (remaining4)
		{
			++cur_num4;
			--remaining4;
		}

		trans_assign.push_back(cur_num5);
		cur_num5 += average5;
		if (remaining5)
		{
			++cur_num5;
			--remaining5;
		}

		dp_assign.push_back(cur_num6);
		cur_num6 += average6;
		if (remaining6)
		{
			++cur_num6;
			--remaining6;
		}
	}

}

void SED_CM_FMM::fillSEDZ_MP()
{
	auto task_num = Task_assign1.size() - 1;
	for (int t = 0; t < task_num; t++)
	{
		int begin = Task_assign1[t];
		int end = Task_assign1[t + 1];
		if (begin < end)
		{
			auto task = std::bind(&SED_CM_FMM::SEDZ_MP_Kernel, this, begin, end);
			pool_.submit(task);
		}
	}
	pool_.interactiveRun();
}

void SED_CM_FMM::SEDZ_MP_Kernel(int begin, int end)
{
	for (int t = begin; t < end; t++)
	{
		auto f = SED_task1[t].first;
		auto s = SED_task1[t].second;

		int start_tf, end_tf, start_ts, end_ts;
		int start_ef, end_ef, start_es, end_es;

		ct_ptr_->getSEDCommonTriangleSize(f, start_tf, end_tf);
		ce_ptr_->getSEDCommonEdgeSize(f, start_ef, end_ef);
		ct_ptr_->getSEDCommonTriangleSize(s, start_ts, end_ts);
		ce_ptr_->getSEDCommonEdgeSize(s, start_es, end_es);



		fillSEDZDD_MP(f, s, start_tf, end_tf, start_ts, end_ts, start_ef, end_ef, start_es, end_es);
		fillSEDZDM(f, s, start_tf, end_tf, start_ts, end_ts, start_ef, end_ef, start_es, end_es);
		fillSEDZMD(f, s, start_tf, end_tf, start_ts, end_ts, start_ef, end_ef, start_es, end_es);
		fillSEDZMM(f, s, start_tf, end_tf, start_ts, end_ts, start_ef, end_ef, start_es, end_es);
	}
}

void SED_CM_FMM::fillREDZ_MP()
{
	auto task_num = Task_assign2.size() - 1;
	Qcout << task_num << std::endl;
	Z_red_eigen.resize(Array_x*Array_y*svd_k, Array_x*Array_y*svd_k);
	prepareSolve();
	for (int t = 0; t < task_num; t++)
	{
		int begin = Task_assign2[t];
		int end = Task_assign2[t + 1];

		//Qcout <<setw(8)<< t << setw(8) << begin << setw(8) << end << std::endl;

		if (begin < end)
		{
			auto task = std::bind(&SED_CM_FMM::REDZ_MP_Kernel, this, begin, end);
			pool_.submit(task);
		}
	}
	Qcout << "Finish task submit!" << std::endl;
	pool_.interactiveRun();
}

void SED_CM_FMM::REDZ_MP_Kernel(int begin, int end)
{
	for (int t = begin; t < end; t++)
	{
		int f = SED_task2[t].first;//x in (x,y)
		int s = SED_task2[t].second;//y in (x,y)

		int x_m = f;
		int y_m = s;

		int x_n = 0;
		int y_n = 0;
		//calculate the transfer factor
		CVector transfer;
		if ((x_m!= 0)||(y_m!=0))
		{
			preCalculateTrans(x_m, y_m, transfer);
		}
		//Qcout << setw(5) << f << setw(5) << s << setw(5) << x_m << setw(5) << y_m << setw(5) << x_n << setw(5) << y_n << setw(5) << id_m << setw(5) << id_n << std::endl;
		//CMatrix Z_coup_all = coupZ_MP(x_m, y_m, x_n, y_n);
		for (int x_offset = 0; x_offset < Array_x; x_offset++)
		{
			int x_temp = f + x_offset;//终点x坐标

			if (x_temp<0)
				continue;
			if (x_temp > Array_x - 1)
				break;

			for (int y_offset = 0; y_offset < Array_y; y_offset++)
			{
				int y_temp = s + y_offset;//终点y坐标

				if (y_temp<0)
					continue;
				if (y_temp > Array_y - 1)
					break;

				int id_m = ArrayToSED(x_temp, y_temp);
				int id_n = ArrayToSED(x_offset, y_offset);

				int start_tf, end_tf, start_ts, end_ts;
				int start_ef, end_ef, start_es, end_es;

				int b_ce1_f, b_ce2_f;
				int b_ce1_s, b_ce2_s;

				ct_ptr_->getSEDCommonTriangleSize(id_m, start_tf, end_tf);
				ce_ptr_->getSEDCommonEdgeSize(id_m, start_ef, end_ef);
				ct_ptr_->getSEDCommonTriangleSize(id_n, start_ts, end_ts);
				ce_ptr_->getSEDCommonEdgeSize(id_n, start_es, end_es);

				getBoundaryNum(id_m, b_ce1_f, b_ce2_f);
				getBoundaryNum(id_n, b_ce1_s, b_ce2_s);
				
				CMatrix Z_coup;
				/*Z_coup.resize(end_ef - start_ef, end_es - start_es);

				Z_coup.block(0, 0, com_ce, com_ce) = Z_coup_all.block(com_ct, com_ct, com_ce, com_ce);

				Z_coup.block(com_ce, 0, b_ce1_f, com_ce) = Z_coup_all.block(com_ce, 0, b_ce1_f, com_ce);
				Z_coup.block(com_ce + b_ce1_f, 0, b_ce2_f, com_ce) = Z_coup_all.block(com_ce + bou_ce1, 0, b_ce2_f, com_ce);

				Z_coup.block(0, com_ce, com_ce, b_ce1_s) = Z_coup_all.block(0, com_ce, com_ce, b_ce1_s);
				Z_coup.block(0, com_ce + b_ce1_s, com_ce, b_ce2_s) = Z_coup_all.block(0, com_ce + bou_ce1, com_ce, b_ce2_s);

				Z_coup.block(com_ce, com_ce, b_ce1_f, b_ce1_s) = Z_coup_all.block(com_ce, com_ce, b_ce1_f, b_ce1_s);
				Z_coup.block(com_ce + b_ce1_f, com_ce, b_ce2_f, b_ce1_s) = Z_coup_all.block(com_ce + bou_ce1, com_ce, b_ce2_f, b_ce1_s);
				Z_coup.block(com_ce, com_ce + b_ce1_s, b_ce1_f, b_ce2_s) = Z_coup_all.block(com_ce, com_ce + bou_ce1, b_ce1_f, b_ce2_s);
				Z_coup.block(com_ce + b_ce1_f, com_ce + b_ce1_s, b_ce2_f, b_ce2_s) = Z_coup_all.block(com_ce + bou_ce1, com_ce + bou_ce1, b_ce2_f, b_ce2_s);*/
				if (x_m == 0 && y_m == 0)//是否是近场自耦合
				{
					Z_coup = coupZ_MP(x_offset, y_offset, x_offset, y_offset);
				}
				for (int k_m = 0; k_m < svd_k; k_m++)
				{
					CVector I_m = I_svd_eigen[id_m].col(k_m);
					int r = (x_temp*Array_y + y_temp) * svd_k + k_m;

					for (int k_n = 0; k_n < svd_k; k_n++)
					{
						CVector I_n = I_svd_eigen[id_n].col(k_n);
						int c = (x_offset*Array_y + y_offset) * svd_k + k_n;

						if (x_m == 0 && y_m == 0)
						{
							CMatrix z_ele = I_m.transpose()*Z_coup*I_n;
							Z_red_eigen.block(r, c, z_ele.rows(), z_ele.cols()) = z_ele;
						}
						else
						{
							CVector V_n = matrixVectorMultiply(x_offset, I_n, transfer);
							CMatrix z_ele = I_m.transpose()*V_n;
							Z_red_eigen.block(r, c, z_ele.rows(), z_ele.cols()) = z_ele;
						}
					}
				}
			}
		}



	}
}

Complex SED_CM_FMM::DDZppKernel(VectorR3 * vf4, VectorR3 * vs4, VectorR3 & vf, VectorR3 & vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_p, value_t &epsp, value_t &epsm)
{
	value_t r, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	//四面体不重合，所以Z1=0

	//z2,z5,z6
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z2,z5
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = exp(-J0 * k_ * r) / r;
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s  *  G);//a2
			Z5 += wt4[p] * wt4[q] * G;
		}

		//z6
		if (abs(epsp - epsm) > 1e-4)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg4f[p] - vg3s[q]).Norm();
				G = exp(-J0 * k_*r) / r;
				Z6 += wt4[p] * w3[q] * G;
			}
		}
	}

	//z3,z4
	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 4; q++)
			{
				r = (vg3f[p] - vg4s[q]).Norm();
				G = exp(-J0 * k_*r) / r;
				Z3 += w3[p] * wt4[q] * G;
			}

			if (abs((epsm - epsp)) > 1e-4)
			{
				for (int q = 0; q < 3; q++)
				{
					r = (vg3f[p] - vg3s[q]).Norm();
					G = exp(-J0 * k_*r) / r;
					Z4 += w3[p] * w3[q] * G;
				}
			}
		}
	}


	Z = (-1.0f*k_* kp / 9.0f) * Z2 + (1.0f / k_) * (-kp * Z3 - (km - kp)*Z4 + kp * Z5 + (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::DDZpmKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_m, value_t &epsp, value_t &epsm)
{
	value_t r, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	//四面体不重合，所以Z1=0

	//z2,z5
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = exp(-J0 * k_ * r) / r;
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s  *  G);//a2
			Z5 += wt4[p] * wt4[q] * G;//a3
		}
	}

	//z3
	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 4; q++)
			{
				r = (vg3f[p] - vg4s[q]).Norm();
				G = exp(-J0 * k_*r) / r;
				Z3 += w3[p] * wt4[q] * G;
			}
		}
	}


	Z = (k_*km / 9.0f)*Z2 + (1.0f / k_) * (km* Z3 - km * Z5);
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::DDZmpKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_p, value_t &epsp, value_t &epsm)
{
	value_t r, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	//四面体不重合，所以Z1=0

	//z2,z5,z6
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = exp(-J0 * k_ * r) / r;
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s  *  G);//a2
			Z5 += wt4[p] * wt4[q] * G;
		}

		if (abs((epsm - epsp)) > 1e-4)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg4f[p] - vg3s[q]).Norm();
				G = exp(-J0 * k_*r) / r;
				Z6 += wt4[p] * w3[q] * G;
			}
		}
	}

	Z = (k_*kp / 9.0f)*Z2 + (1.0f / k_) *(-kp * Z5 - (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::DDZmmKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_m, value_t &epsp, value_t &epsm)
{
	value_t r, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	//四面体不重合，所以Z1=0

	//z2,z5
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = exp(-J0 * k_ * r) / r;
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s  *  G);//a2
			Z5 += wt4[p] * wt4[q] * G;
		}
	}



	Z = (-k_ * km / 9.0f)*Z2 + (1.0f / k_) *(km*Z5);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::DDZppSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_p, value_t &epsp, value_t &epsm)
{
	value_t r, volume, area, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4s[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}
	//z2,z5,z6
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z2,z5非奇异部分
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s  *  G);
			Z5 += wt4[p] * wt4[q] * G;
		}
		//z2,z5奇异部分
		VectorR3 Ivp(0.0f, 0.0f, 0.0f);
		value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
			Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
			Iv += 0.5f*IvSingular(v4_3[j], vg4f[p]);
		}
		Z2 += wt4[p] * (rou_f ^ (Ivp + Iv * (vg4f[p] - vs))) / volume;
		Z5 += wt4[p] * (Iv / volume);

		//z6
		if (abs(epsp - epsm) > 1e-4)
		{
			//z6非奇异部分
			for (int q = 0; q < 3; q++)
			{
				r = (vg4f[p] - vg3s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z6 += wt4[p] * w3[q] * G;
			}
			//z6奇异部分
			Z6 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
		}
	}
	//z3,z4非奇异部分
	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 4; q++)
			{
				r = (vg3f[p] - vg4s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z3 += w3[p] * wt4[q] * G;
			}
			//z3奇异部分
			value_t Iv = 0.0f;
			for (int j = 0; j < 4; j++)
			{
				Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
			}
			Z3 += w3[p] * Iv / volume;

			//z4非奇异
			if (abs(epsp - epsm) > 1e-4)
			{
				for (int q = 0; q < 3; q++)
				{
					r = (vg3f[p] - vg3s[q]).Norm();
					G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
					Z4 += w3[p] * w3[q] * G;
				}
				//z4奇异
				Z4 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
			}
		}
	}
	Z = (1.0f / k_) *(PI4 / (9.0f*volume*epsp))*Z1 - (k_* kp / 9.0f)*Z2 + (1.0f / k_) *(-kp * Z3 - (km - kp)*Z4 + kp * Z5 + (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::DDZpmSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_m, value_t &epsp, value_t &epsm)
{
	value_t r, volume, area, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4s[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}
	//z2,z5
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z2,z5非奇异部分
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s  *  G);
			Z5 += wt4[p] * wt4[q] * G;
		}
		//z2,z5奇异部分
		VectorR3 Ivp(0, 0, 0); value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
			Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
			Iv += 0.5f*IvSingular(v4_3[j], vg4f[p]);
		}
		Z2 += wt4[p] * rou_f ^ (Ivp + Iv * (vg4f[p] - vs)) / volume;
		Z5 += wt4[p] * Iv / volume;

		//z6
	}
	//z3,z4非奇异部分
	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 4; q++)
			{
				r = (vg3f[p] - vg4s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z3 += w3[p] * wt4[q] * G;
			}
			//z3奇异部分
			value_t Iv = 0.0f;
			for (int j = 0; j < 4; j++)
			{
				Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
			}
			Z3 += w3[p] * Iv / volume;

			//z4
		}
	}

	Z = (-1.0f / k_) *(PI4 / (9.0f * volume * epsm))*Z1 + (k_*km / 9.0f)*Z2 + (1.0f / k_) *(km*Z3 - km * Z5);
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::DDZmpSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_p, value_t &epsp, value_t &epsm)
{
	value_t r, volume, area, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//四面体不重合，所以Z1=0
	//z1
	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4s[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}
	//z2,z5,z6
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z2,z5非奇异部分
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s  *  G);
			Z5 += wt4[p] * wt4[q] * G;
		}
		//z2,z5奇异部分
		VectorR3 Ivp(0, 0, 0); value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
			Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
			Iv += 0.5f*IvSingular(v4_3[j], vg4f[p]);
		}
		Z2 += wt4[p] * rou_f ^ (Ivp + Iv * (vg4f[p] - vs)) / volume;
		Z5 += wt4[p] * Iv / volume;


		if (abs(epsp - epsm) > 1e-4)
		{
			//z6非奇异部分
			for (int q = 0; q < 3; q++)
			{
				r = (vg4f[p] - vg3s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z6 += wt4[p] * w3[q] * G;
			}
			//z6奇异部分
			value_t Is = 0.0f;
			Z6 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
		}
	}
	//z3,z4非奇异部分

	Z = (-1.0f / k_) *(PI4 / (9.0f*volume*epsp))*Z1 + (k_*kp / 9.0f)*Z2 + (1.0f / k_) *(-kp * Z5 - (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::DDZmmSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_m, value_t &epsp, value_t &epsm)
{
	value_t r, volume, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//四面体不重合，所以Z1=0
	//z1
	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4s[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}
	//z2,z5,z6
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z2,z5非奇异部分
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s  *  G);
			Z5 += wt4[p] * wt4[q] * G;
		}
		//z2,z5奇异部分
		VectorR3 Ivp(0, 0, 0); value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
			Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
			Iv += 0.5f*IvSingular(v4_3[j], vg4f[p]);
		}
		Z2 += wt4[p] * rou_f ^ (Ivp + Iv * (vg4f[p] - vs)) / volume;
		Z5 += wt4[p] * Iv / volume;
	}
	//z3,z4非奇异部分


	Z = (1.0f / k_) *(PI4 / (9.0f*volume*epsm))*Z1 - (k_* km / 9.0f)*Z2 + (1.0f / k_) *(km*Z5);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}


Complex SED_CM_FMM::DDZppKernel_fast(VectorR3 * vf4, VectorR3 * vs4, VectorR3 & vf, VectorR3 & vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_p, value_t &epsp, value_t &epsm, KIM& ik)
{
	value_t r, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	//四面体不重合，所以Z1=0

	//z2,z5
	Z2 = ik.Z2K1(t_f_p, t_s_p) - (vf^ik.Z2K2(t_f_p, t_s_p)) - (vs^ik.Z2K3(t_f_p, t_s_p)) + ((vs^vf)*ik.Z2K4(t_f_p, t_s_p));
	Z5 = ik.Z5K(t_f_p, t_s_p);

	//z6
	for (int p = 0; p < 4; p++)
	{
		if (abs(epsp - epsm) > 1e-4)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg4f[p] - vg3s[q]).Norm();
				G = exp(-J0 * k_*r) / r;
				Z6 += wt4[p] * w3[q] * G;
			}
		}
	}

	//z3,z4
	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 4; q++)
			{
				r = (vg3f[p] - vg4s[q]).Norm();
				G = exp(-J0 * k_*r) / r;
				Z3 += w3[p] * wt4[q] * G;
			}

			if (abs((epsm - epsp)) > 1e-4)
			{
				for (int q = 0; q < 3; q++)
				{
					r = (vg3f[p] - vg3s[q]).Norm();
					G = exp(-J0 * k_*r) / r;
					Z4 += w3[p] * w3[q] * G;
				}
			}
		}
	}


	Z = (-1.0f*k_* kp / 9.0f) * Z2 + (1.0f / k_) * (-kp * Z3 - (km - kp)*Z4 + kp * Z5 + (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::DDZpmKernel_fast(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_m, value_t &epsp, value_t &epsm, KIM& ik)
{
	value_t r, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	//四面体不重合，所以Z1=0

	//z2,z5
	Z2 = ik.Z2K1(t_f_p, t_s_m) - (vf^ik.Z2K2(t_f_p, t_s_m)) - (vs^ik.Z2K3(t_f_p, t_s_m)) + ((vs^vf)*ik.Z2K4(t_f_p, t_s_m));
	Z5 = ik.Z5K(t_f_p, t_s_m);

	//z3
	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 4; q++)
			{
				r = (vg3f[p] - vg4s[q]).Norm();
				G = exp(-J0 * k_*r) / r;
				Z3 += w3[p] * wt4[q] * G;
			}
		}
	}


	Z = (k_*km / 9.0f)*Z2 + (1.0f / k_) * (km* Z3 - km * Z5);
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::DDZmpKernel_fast(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_p, value_t &epsp, value_t &epsm, KIM& ik)
{
	value_t r, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	//四面体不重合，所以Z1=0

	//z2,z5
	Z2 = ik.Z2K1(t_f_m, t_s_p) - (vf^ik.Z2K2(t_f_m, t_s_p)) - (vs^ik.Z2K3(t_f_m, t_s_p)) + ((vs^vf)*ik.Z2K4(t_f_m, t_s_p));
	Z5 = ik.Z5K(t_f_m, t_s_p);

	//z6
	for (int p = 0; p < 4; p++)
	{
		if (abs((epsm - epsp)) > 1e-4)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg4f[p] - vg3s[q]).Norm();
				G = exp(-J0 * k_*r) / r;
				Z6 += wt4[p] * w3[q] * G;
			}
		}
	}

	Z = (k_*kp / 9.0f)*Z2 + (1.0f / k_) *(-kp * Z5 - (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::DDZmmKernel_fast(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_m, value_t &epsp, value_t &epsm, KIM& ik)
{
	value_t r, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	//四面体不重合，所以Z1=0

	//z2,z5
	Z2 = ik.Z2K1(t_f_m, t_s_m) - (vf^ik.Z2K2(t_f_m, t_s_m)) - (vs^ik.Z2K3(t_f_m, t_s_m)) + ((vs^vf)*ik.Z2K4(t_f_m, t_s_m));
	Z5 = ik.Z5K(t_f_m, t_s_m);



	Z = (-k_ * km / 9.0f)*Z2 + (1.0f / k_) *(km*Z5);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::DDZppSingular_fast(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_p, value_t &epsp, value_t &epsm, KIM& ik)
{
	value_t r, volume, area, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//四面体不重合，所以Z1=0
	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4s[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}
	//z2,z5,z6
	Z2 = ik.Z2K1(t_f_p, t_s_p) - (vf^ik.Z2K2(t_f_p, t_s_p)) - (vs^ik.Z2K3(t_f_p, t_s_p)) + ((vs^vf)*ik.Z2K4(t_f_p, t_s_p)) +
		((ik.RIVP(t_f_p) + ik.RRIV(t_f_p) - (vs^ik.RIV[t_f_p]) - (vf^ik.IVP[t_f_p]) - (vf^ik.RIV[t_f_p]) + ((vf^vs)*ik.IV(t_f_p))) / volume);
	Z5 = ik.Z5K(t_f_p, t_s_p) + (ik.IV(t_f_p) / volume);

	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z6
		if (abs(epsp - epsm) > 1e-4)
		{
			//z6非奇异部分
			for (int q = 0; q < 3; q++)
			{
				r = (vg4f[p] - vg3s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z6 += wt4[p] * w3[q] * G;
			}
			//z6奇异部分
			Z6 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
		}
	}
	//z3,z4非奇异部分
	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 4; q++)
			{
				r = (vg3f[p] - vg4s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z3 += w3[p] * wt4[q] * G;
			}
			//z3奇异部分
			value_t Iv = 0.0f;
			for (int j = 0; j < 4; j++)
			{
				Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
			}
			Z3 += w3[p] * Iv / volume;

			//z4非奇异
			if (abs(epsp - epsm) > 1e-4)
			{
				for (int q = 0; q < 3; q++)
				{
					r = (vg3f[p] - vg3s[q]).Norm();
					G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
					Z4 += w3[p] * w3[q] * G;
				}
				//z4奇异
				Z4 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
			}
		}
	}
	Z = (1.0f / k_) *(PI4 / (9.0f*volume*epsp))*Z1 - (k_* kp / 9.0f)*Z2 + (1.0f / k_) *(-kp * Z3 - (km - kp)*Z4 + kp * Z5 + (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::DDZpmSingular_fast(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_m, value_t &epsp, value_t &epsm, KIM& ik)
{
	value_t r, volume, area, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4s[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}
	//z2,z5
	Z2 = ik.Z2K1(t_f_p, t_s_m) - (vf^ik.Z2K2(t_f_p, t_s_m)) - (vs^ik.Z2K3(t_f_p, t_s_m)) + ((vs^vf)*ik.Z2K4(t_f_p, t_s_m)) +
		((ik.RIVP(t_f_p) + ik.RRIV(t_f_p) - (vs^ik.RIV[t_f_p]) - (vf^ik.IVP[t_f_p]) - (vf^ik.RIV[t_f_p]) + ((vf^vs)*ik.IV(t_f_p))) / volume);
	Z5 = ik.Z5K(t_f_p, t_s_m) + (ik.IV(t_f_p) / volume);

	//z3,z4非奇异部分
	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 4; q++)
			{
				r = (vg3f[p] - vg4s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z3 += w3[p] * wt4[q] * G;
			}
			//z3奇异部分
			value_t Iv = 0.0f;
			for (int j = 0; j < 4; j++)
			{
				Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
			}
			Z3 += w3[p] * Iv / volume;

			//z4
		}
	}

	Z = (-1.0f / k_) *(PI4 / (9.0f * volume * epsm))*Z1 + (k_*km / 9.0f)*Z2 + (1.0f / k_) *(km*Z3 - km * Z5);
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::DDZmpSingular_fast(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_p, value_t &epsp, value_t &epsm, KIM& ik)
{
	value_t r, volume, area, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//四面体不重合，所以Z1=0
	//z1
	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4s[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}
	//z2,z5
	Z2 = ik.Z2K1(t_f_m, t_s_p) - (vf^ik.Z2K2(t_f_m, t_s_p)) - (vs^ik.Z2K3(t_f_m, t_s_p)) + ((vs^vf)*ik.Z2K4(t_f_m, t_s_p)) +
		((ik.RIVP(t_f_m) + ik.RRIV(t_f_m) - (vs^ik.RIV[t_f_m]) - (vf^ik.IVP[t_f_m]) - (vf^ik.RIV[t_f_m]) + ((vf^vs)*ik.IV(t_f_m))) / volume);
	Z5 = ik.Z5K(t_f_m, t_s_p) + (ik.IV(t_f_m) / volume);

	//z6
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		if (abs(epsp - epsm) > 1e-4)
		{
			//z6非奇异部分
			for (int q = 0; q < 3; q++)
			{
				r = (vg4f[p] - vg3s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z6 += wt4[p] * w3[q] * G;
			}
			//z6奇异部分
			value_t Is = 0.0f;
			Z6 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
		}
	}
	//z3,z4非奇异部分

	Z = (-1.0f / k_) *(PI4 / (9.0f*volume*epsp))*Z1 + (k_*kp / 9.0f)*Z2 + (1.0f / k_) *(-kp * Z5 - (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::DDZmmSingular_fast(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_m, value_t &epsp, value_t &epsm, KIM& ik)
{
	value_t r, volume, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//四面体不重合，所以Z1=0
	//z1
	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4s[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}
	//z2,z5
	Z2 = ik.Z2K1(t_f_m, t_s_m) - (vf^ik.Z2K2(t_f_m, t_s_m)) - (vs^ik.Z2K3(t_f_m, t_s_m)) + ((vs^vf)*ik.Z2K4(t_f_m, t_s_m)) +
		((ik.RIVP(t_f_m) + ik.RRIV(t_f_m) - (vs^ik.RIV[t_f_m]) - (vf^ik.IVP[t_f_m]) - (vf^ik.RIV[t_f_m]) + ((vf^vs)*ik.IV(t_f_m))) / volume);
	Z5 = ik.Z5K(t_f_m, t_s_m) + (ik.IV(t_f_m) / volume);


	Z = (1.0f / k_) *(PI4 / (9.0f*volume*epsm))*Z1 - (k_* km / 9.0f)*Z2 + (1.0f / k_) *(km*Z5);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}


Complex SED_CM_FMM::DMZppKernel(VectorR3 * vf4, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 * vf3, int & _vm_f)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg3f[3], vg3s[3], rou_f, rou_s;
	value_t r;
	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);

	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vg3s[q] - vs;
			r = (vg4f[p] - vg3s[q]).Norm();
			G = exp(-J0 * k_*r) / r;
			Z1 += wt4[p] * w3[q] * (rou_f^rou_s)*G;
			Z2 += wt4[p] * w3[q] * G;
		}
	}

	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg3f[p] - vg3s[q]).Norm();
				G = exp(-J0 * k_*r) / r;
				Z3 += w3[p] * w3[q] * G;
			}
		}
	}

	Z = (1.0f / 6.0f)*Z1 - (1.0f / (k_*k_))*(Z2 - Z3);//外面乘以J0*k_*Z0/PI4

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}

	return Z;
}

Complex SED_CM_FMM::DMZpmKernel(VectorR3 * vf4, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 * vf3, int & _vm_f)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg3f[3], vg3s[3], rou_f, rou_s;
	value_t r;
	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);

	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vg3s[q] - vs;
			r = (vg4f[p] - vg3s[q]).Norm();
			G = exp(-J0 * k_*r) / r;
			Z1 += wt4[p] * w3[q] * (rou_f^rou_s)*G;
			Z2 += wt4[p] * w3[q] * G;
		}
	}

	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg3f[p] - vg3s[q]).Norm();
				G = exp(-J0 * k_*r) / r;
				Z3 += w3[p] * w3[q] * G;
			}
		}
	}


	Z = -(1.0f / 6.0f)*Z1 + (1.0f / (k_*k_))*(Z2 - Z3);//外面乘以J0*k_*Z0

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::DMZmpKernel(VectorR3 * vf4, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 * vf3, int & _vm_f)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg3f[3], vg3s[3], rou_f, rou_s;
	value_t r;
	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);

	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vg3s[q] - vs;
			r = (vg4f[p] - vg3s[q]).Norm();
			G = exp(-J0 * k_*r) / r;
			Z1 += wt4[p] * w3[q] * (rou_f^rou_s)*G;
			Z2 += wt4[p] * w3[q] * G;
		}
	}

	/*if (_vm_f == -1)
	{
	for (int p = 0; p < 3; p++)
	{
	for (int q = 0; q < 3; q++)
	{
	r = (vg1f[p] - vg3s[q]).Norm();
	G = exp(-J0*k_*r) / r;
	Z3 += w3[p] * w3[q] * G;
	}
	}
	}*/

	Z = -(1.0f / 6.0f)*Z1 + (1.0f / (k_*k_))*(Z2);//外面乘以J0*k_*Z0

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::DMZmmKernel(VectorR3 * vf4, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 * vf3, int & _vm_f)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg3f[3], vg3s[3], rou_f, rou_s;
	value_t r;
	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);

	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vg3s[q] - vs;
			r = (vg4f[p] - vg3s[q]).Norm();
			G = exp(-J0 * k_*r) / r;
			Z1 += wt4[p] * w3[q] * (rou_f^rou_s)*G;
			Z2 += wt4[p] * w3[q] * G;
		}
	}

	/*if (_vm_f == -1)
	{
	for (int p = 0; p < 3; p++)
	{
	for (int q = 0; q < 3; q++)
	{
	r = (vg1f[p] - vg3s[q]).Norm();
	G = exp(-J0*k_*r) / r;
	Z3 += w3[p] * w3[q] * G;
	}
	}
	}*/

	Z = (1.0f / 6.0f)*Z1 - (1.0f / (k_*k_))*(Z2);//外面乘以J0*k_*Z0

	return Z;
}


Complex	SED_CM_FMM::DMZppSingular(VectorR3 * vf4, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 * vf3, int & _vm_f)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg3f[3], vg3s[3], rou_f, rou_s;
	value_t r, area;
	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	area = Area(vs3[0], vs3[1], vs3[2]);

	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//非奇异部分
		for (int q = 0; q < 3; q++)
		{
			rou_s = vg3s[q] - vs;
			r = (vg4f[p] - vg3s[q]).Norm();
			G = exp(-J0 * k_*r) / r;
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z1 += wt4[p] * w3[q] * (rou_f^rou_s)*G;
			Z2 += wt4[p] * w3[q] * G;
		}
		//奇异部分
		//Z1 += wt4[p] * (rou_f ^ (IspSingular(vs3, vg4f[p]) + (vg4f[p] - vs)*IsSingular(vs3, vg4f[p]))) / area;
		//Z2 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
	}

	for (int p = 0; p < 3; p++)
	{
		//Z3非奇异部分
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z3 += w3[p] * w3[q] * G;
		}
		//Z3奇异部分
		Z3 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	Z = (1.0f / 6.0f)*Z1 - (1.0f / (k_*k_))*(Z2 - Z3);//外面乘以J0*k_*Z0/PI4

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}

	return Z;
}

Complex	SED_CM_FMM::DMZpmSingular(VectorR3 * vf4, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 * vf3, int & _vm_f)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg3f[3], vg3s[3], rou_f, rou_s;
	value_t r, area;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	area = Area(vs3[0], vs3[1], vs3[2]);

	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//非奇异部分
		for (int q = 0; q < 3; q++)
		{
			rou_s = vg3s[q] - vs;
			r = (vg4f[p] - vg3s[q]).Norm();
			G = exp(-J0 * k_*r) / r;
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z1 += wt4[p] * w3[q] * (rou_f^rou_s)*G;
			Z2 += wt4[p] * w3[q] * G;
		}
		//奇异部分
		//Z1 += wt4[p] * (rou_f ^ (IspSingular(vs3, vg4f[p]) + (vg4f[p] - vs)*IsSingular(vs3, vg4f[p]))) / area;
		//Z2 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
	}

	for (int p = 0; p < 3; p++)
	{
		//Z3非奇异部分
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z3 += w3[p] * w3[q] * G;
		}
		//Z3奇异部分
		Z3 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	Z = -(1.0f / 6.0f)*Z1 + (1.0f / (k_*k_))*(Z2 - Z3);//外面乘以J0*k_*Z0/PI4

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}

	return Z;
}


Complex SED_CM_FMM::MDZppkernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg3f[3], vg3s[3], vg4s[4], rou_f, rou_s;
	value_t r, kp, km;
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	kp = 1.0f - (1.0f / epsp);
	km = 1.0f - (1.0f / epsm);
	for (int p = 0; p < 3; p++)
	{
		rou_f = vg3f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg3f[p] - vg4s[q]).Norm();
			G = exp(-J0 * k_*r) / r;
			Z1 += w3[p] * wt4[q] * (rou_f^rou_s)*G;
			Z2 += w3[p] * wt4[q] * G;//in the end need *kp
		}
	}

	if (abs(epsm - epsp)>1e-6)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg3f[p] - vg3s[q]).Norm();
				G = exp(-J0 * k_*r) / r;
				Z3 += w3[p] * w3[q] * G;//in the end need *(km-kp)
			}
		}
	}

	Z = (1.0f / 6.0f)*(-kp * Z1 + (6.0f / (k_*k_))*(kp*Z2 + (km - kp)*Z3));

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::MDZpmkernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg3f[3], vg3s[3], vg4s[4], rou_f, rou_s;
	value_t r, kp, km;
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	kp = 1.0f - (1.0f / epsp);
	km = 1.0f - (1.0f / epsm);
	for (int p = 0; p < 3; p++)
	{
		rou_f = vg3f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg3f[p] - vg4s[q]).Norm();
			G = exp(-J0 * k_*r) / r;
			Z1 += w3[p] * wt4[q] * (rou_f^rou_s)*G;
			Z2 += w3[p] * wt4[q] * G;//in the end need *kp
		}
	}

	/*if (abs(epsm - epsp)<1e-6)
	{
	for (int p = 0; p < 3; p++)
	{
	for (int q = 0; q < 3; q++)
	{
	r = (vg3f[p] - vg5s[q]).Norm();
	G = exp(-J0*k_*r) / r;
	Z3 += w3[p] * w3[q] * G;//in the end need *(km-kp)
	}
	}
	}*/

	Z = (1.0f / 6.0f)*(km*Z1 - (6.0f / (k_*k_))*(km*Z2 + (km - kp)*Z3));
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::MDZmpkernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg3f[3], vg3s[3], vg4s[4], rou_f, rou_s;
	value_t r, kp, km;
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	kp = 1.0f - (1.0f / epsp);
	km = 1.0f - (1.0f / epsm);
	for (int p = 0; p < 3; p++)
	{
		rou_f = vg3f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg3f[p] - vg4s[q]).Norm();
			G = exp(-J0 * k_*r) / r;
			Z1 += w3[p] * wt4[q] * (rou_f^rou_s)*G;
			Z2 += w3[p] * wt4[q] * G;//in the end need *kp
		}
	}

	if (abs(epsm - epsp)>1e-6)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg3f[p] - vg3s[q]).Norm();
				G = exp(-J0 * k_*r) / r;
				Z3 += w3[p] * w3[q] * G;//in the end need *(km-kp)
			}
		}
	}

	Z = (1.0f / 6.0f)*(kp*Z1 - (6.0f / (k_*k_))*(kp*Z2 + (km - kp)*Z3));
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::MDZmmkernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg3f[3], vg3s[3], vg4s[4], rou_f, rou_s;
	value_t r, kp, km;
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	kp = 1.0f - (1.0f / epsp);
	km = 1.0f - (1.0f / epsm);
	for (int p = 0; p < 3; p++)
	{
		rou_f = vg3f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg3f[p] - vg4s[q]).Norm();
			G = exp(-J0 * k_*r) / r;
			Z1 += w3[p] * wt4[q] * (rou_f^rou_s)*G;
			Z2 += w3[p] * wt4[q] * G;//in the end need *kp
		}
	}

	/*if (abs(epsm - epsp)<1e-6)
	{
	for (int p = 0; p < 3; p++)
	{
	for (int q = 0; q < 3; q++)
	{
	r = (vg3f[p] - vg5s[q]).Norm();
	G = exp(-J0*k_*r) / r;
	Z3 += w3[p] * w3[q] * G;//in the end need *(km-kp)
	}
	}
	}*/


	Z = (1.0f / 6.0f)*(-km * Z1 + (6.0f / (k_*k_))*(km*Z2 + (km - kp)*Z3));
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}


Complex SED_CM_FMM::MDZppSingular(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg3f[3], vg3s[3], vg4s[4], rou_f, rou_s;
	value_t r, kp, km;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	kp = 1.0f - (1.0f / epsp);
	km = 1.0f - (1.0f / epsm);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);
	value_t volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vg3f[p] - vf;
		//非奇异部分
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg3f[p] - vg4s[q]).Norm();
			G = exp(-J0 * k_*r) / r;
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z1 += w3[p] * wt4[q] * (rou_f^rou_s)*G;
			Z2 += w3[p] * wt4[q] * G;//in the end need *kp
		}
		//奇异部分
		/*VectorR3 Ivp(0.0f, 0.0f, 0.0f);
		value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
		Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg3f[p]);
		Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
		}
		Z1 += w3[p] * (rou_f ^ (Ivp + Iv*(vg3f[p] - vs))) / volume;
		Z2 += w3[p] * (Iv / volume);*/
	}

	value_t area = Area(vs3[0], vs3[1], vs3[2]);
	for (int p = 0; p < 3; p++)
	{
		//Z3非奇异部分
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z3 += w3[p] * w3[q] * G;//in the end need *(km-kp)
		}
		//Z3奇异部分
		Z3 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	Z = (1.0f / 6.0f)*(-kp * Z1 + (6.0f / (k_*k_))*(kp*Z2 + (km - kp)*Z3));
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex SED_CM_FMM::MDZmpSingular(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg3f[3], vg3s[3], vg4s[4], rou_f, rou_s;
	value_t r, kp, km;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	kp = 1.0f - (1.0f / epsp);
	km = 1.0f - (1.0f / epsm);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);
	value_t volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vg3f[p] - vf;
		//非奇异部分
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg3f[p] - vg4s[q]).Norm();
			G = exp(-J0 * k_*r) / r;
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z1 += w3[p] * wt4[q] * (rou_f^rou_s)*G;
			Z2 += w3[p] * wt4[q] * G;//in the end need *kp
		}
		//奇异部分
		/*VectorR3 Ivp(0.0f, 0.0f, 0.0f);
		value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
		Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg3f[p]);
		Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
		}
		Z1 += w3[p] * (rou_f ^ (Ivp + Iv*(vg3f[p] - vs))) / volume;
		Z2 += w3[p] * (Iv / volume);*/
	}

	value_t area = Area(vs3[0], vs3[1], vs3[2]);
	for (int p = 0; p < 3; p++)
	{
		//Z3非奇异部分
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z3 += w3[p] * w3[q] * G;//in the end need *(km-kp)
		}
		//Z3奇异部分
		Z3 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}


	Z = (1.0f / 6.0f)*(kp*Z1 - (6.0f / (k_*k_))*(kp*Z2 + (km - kp)*Z3));
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}


Complex SED_CM_FMM::MMZppKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	value_t r;
	Complex Ze(0, 0), Zm(0, 0), G(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs;
			r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k_ * r) / r;
			Ze += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * G;
		}
	}
	return Ze;
}

Complex SED_CM_FMM::MMZpmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	value_t r;
	Complex Ze(0, 0), Zm(0, 0), G(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vgs[q];
			r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k_ * r) / r;
			Ze += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * G;
		}
	}
	return Ze;
}

Complex SED_CM_FMM::MMZmpKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	value_t r;
	Complex Ze(0, 0), Zm(0, 0), G(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vf - vgf[p];
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs;
			r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k_ * r) / r;
			Ze += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * G;
		}
	}
	return Ze;
}

Complex SED_CM_FMM::MMZmmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	value_t r;
	Complex Ze(0, 0), Zm(0, 0), G(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vf - vgf[p];
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vgs[q];
			r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k_ * r) / r;
			Ze += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * G;
		}
	}
	return Ze;
}

Complex SED_CM_FMM::MMZppSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	/*Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0);// Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	VectorR3 vsc = (vs3[0] + vs3[1] + vs3[2]) / 3.0;
	VectorR3 rou_sc = vsc - vs; //.1.

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
	//电场  非奇异部分
	rou_f = vgf[p] - vf;
	for (int q = 0; q < 3; q++)
	{
	rou_s = vgs[q] - vs; //.2.
	value_t r = (vgf[p] - vgs[q]).Norm();
	Complex coff = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
	Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * coff; //.3.
	}
	//电场  奇异部分
	Complex Ztemp(0, 0);
	value_t Is2(0.0f);
	VectorR3 Is1(0.0f, 0.0f, 0.0f);
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[1]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vgf[p] - ni*(ni^vgf[p]);
	VectorR3 rou_n = vs - ni*(ni^vs);
	for (int i = 0; i < 3; i++)
	{
	VectorR3 li = vs3[(i + 1) % 3] - vs3[i];//not normalize
	VectorR3 lj = vs3[(i + 2) % 3] - vs3[(i + 1) % 3];//not normalize
	VectorR3 r_f = vgf[p] - vs3[i];
	VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
	VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
	VectorR3 u = li*n;
	VectorR3 ui = u.Normalize();//normalize
	VectorR3 li_1 = li.Normalize();//normalize
	value_t p0 = abs((vs3[i] - vgf[p]) ^ ui);
	value_t d = (vgf[p] - vs3[i]) ^ ni;
	//value_t lp = (rou_p - rou) ^ li_1;
	//value_t lm = (rou_m - rou) ^ li_1;
	value_t lp = (vs3[(i + 1) % 3] - vgf[p]) ^ li_1;
	value_t lm = (vs3[i] - vgf[p]) ^ li_1;
	//value_t rp = sqrt(p0*p0 + d*d + lp*lp);
	//value_t rm = sqrt(p0*p0 + d*d + lm*lm);
	value_t rp = (vgf[p] - vs3[(i + 1) % 3]).Norm();
	value_t rm = (vgf[p] - vs3[i]).Norm();

	if (p0 < 1e-10)
	continue;
	VectorR3 pi = ((rou_p - rou) - lp * li_1).Normalize();
	value_t r02 = p0*p0 + d*d;
	Is1 += (0.5f*ui)*(r02*log((rp + lp) / (rm + lm)) + (lp*rp) - (lm*rm));
	Is2 += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
	}
	Ze2 += w3[p] * (0.25f * (rou_f ^ (Is1 + Is2*(rou - rou_n))) - (Is2 / (k_ * k_)));
	}*/

	/*if (abs(Sigma) > 1e-8)
	{
	for (int i = 0; i < 3; i++)
	{
	rou_f = vgf[i] - vf;
	rou_s = vgf[i] - vs;
	Ze3 += w3[i] * rou_f^rou_s;
	}
	Ze3 = (PI / (S_s*J0*k_*Z0*Sigma))*Ze3;
	}*/
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	VectorR3 vsc = (vs3[0] + vs3[1] + vs3[2]) / 3.0;
	VectorR3 rou_sc = vsc - vs; //.1.

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
		//电场  非奇异部分
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs; //.2.
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * coff; //.3.
		}
		//电场  奇异部分
		Complex Ztemp(0, 0);
		for (int i = 0; i < 3; i++)
		{
			VectorR3 Li = vs3[(i + 1) % 3] - vs3[i];
			Li.Normalize();//边向量

			value_t Lpi = (vs3[(i + 1) % 3] - vgf[p]) ^ Li; //(边正顶点 - 场点)·(边向量)
			value_t Lmi = (vs3[i] - vgf[p]) ^ Li; //(边负顶点 - 场点)·(边向量)
			value_t Rpi = (vs3[(i + 1) % 3] - vgf[p]).Norm(); //边正顶点 到 场点 的距离
			value_t Rmi = (vs3[i % 3] - vgf[p]).Norm(); //边负顶点 到 场点 的距离
			VectorR3 sj0 = vs3[i] + Li * ((vgf[p] - vs3[i]) ^ Li); //垂足在u方向位置矢量
			value_t P0i = (sj0 - vgf[p]).Norm(); //垂足 到 场点 的距离
			if (P0i > 1e-10)
				Ztemp += P0i * log((Rpi + Lpi) / (Rmi + Lmi)); //.4.
		}
		Ze2 += w3[p] * (0.25f * (rou_f ^ rou_sc) - 1.0f / (k_ * k_))*Ztemp;
	}

	Ze = (Ze1 + 1.0f / S_s * Ze2);// +Ze3;

	return Ze;
}

Complex SED_CM_FMM::MMZpmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	/*Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0);// Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	VectorR3 vsc = (vs3[0] + vs3[1] + vs3[2]) / 3.0;
	VectorR3 rou_sc = vs - vsc; //.1.

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
	//电场  非奇异部分
	rou_f = vgf[p] - vf;
	for (int q = 0; q < 3; q++)
	{
	rou_s = vs - vgs[q]; //.2.
	value_t r = (vgf[p] - vgs[q]).Norm();
	Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
	Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * coff; //.3.
	}
	//电场  奇异部分

	Complex Ztemp(0.0, 0.0);
	value_t Is2(0.0f);
	VectorR3 Is1(0.0f, 0.0f, 0.0f);
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[1]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vgf[p] - ni*(ni^vgf[p]);
	VectorR3 rou_n = vs - ni*(ni^vs);
	for (int i = 0; i < 3; i++)
	{
	VectorR3 li = vs3[(i + 1) % 3] - vs3[i];//not normalize
	VectorR3 lj = vs3[(i + 2) % 3] - vs3[(i + 1) % 3];//not normalize
	VectorR3 r_f = vgf[p] - vs3[i];
	VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
	VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
	VectorR3 u = li*n;
	VectorR3 ui = u / u.Norm();//normalize
	VectorR3 li_1 = li / li.Norm();//normalize
	value_t p0 = abs((vs3[i] - vgf[p]) ^ ui);
	value_t d = (vgf[p] - vs3[i]) ^ ni;
	value_t lp = (rou_p - rou) ^ li_1;
	value_t lm = (rou_m - rou) ^ li_1;
	value_t rp = sqrt(p0*p0 + d*d + lp*lp);
	value_t rm = sqrt(p0*p0 + d*d + lm*lm);

	if (p0 < 1e-10)
	continue;
	VectorR3 pi = ((rou_p - rou) - lp * li_1).Normalize();
	value_t r02 = p0*p0 + d*d;
	Is1 += (0.5f*ui)*(r02*log((rp + lp) / (rm + lm)) + (lp*rp) - (lm*rm));
	Is2 += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
	}
	Ze2 += w3[p] * (0.25f * (-1.0f*rou_f ^ (Is1 + Is2*(rou - rou_n))) + (Is2 / (k_ * k_)));
	}*/

	/*if (abs(Sigma) > 1e-8)
	{
	for (int i = 0; i < 3; i++)
	{
	rou_f = vgf[i] - vf;
	rou_s = vgf[i] - vs;
	Ze3 += w3[i] * rou_f^rou_s;
	}
	Ze3 = (PI / (S_s*J0*k_*Z0*Sigma))*Ze3;
	}*/
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	VectorR3 vsc = (vs3[0] + vs3[1] + vs3[2]) / 3.0;
	VectorR3 rou_sc = vs - vsc; //.1.

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
		//电场  非奇异部分
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vgs[q]; //.2.
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * coff; //.3.
		}
		//电场  奇异部分

		Complex Ztemp(0.0, 0.0);
		for (int i = 0; i < 3; i++)
		{
			VectorR3 Li = vs3[(i + 1) % 3] - vs3[i];
			Li.Normalize();//边向量

			value_t Lpi = (vs3[(i + 1) % 3] - vgf[p]) ^ Li; //(边正顶点 - 场点)·(边向量)
			value_t Lmi = (vs3[i] - vgf[p]) ^ Li; //(边负顶点 - 场点)·(边向量)
			value_t Rpi = (vs3[(i + 1) % 3] - vgf[p]).Norm(); //边正顶点 到 场点 的距离
			value_t Rmi = (vs3[i % 3] - vgf[p]).Norm(); //边负顶点 到 场点 的距离
			VectorR3 sj0 = vs3[i] + Li * ((vgf[p] - vs3[i]) ^ Li); //垂足在u方向位置矢量
			value_t P0i = (sj0 - vgf[p]).Norm(); //垂足 到 场点 的距离
			if (P0i > 1e-10)
				Ztemp += P0i * log((Rpi + Lpi) / (Rmi + Lmi)); //.4.
		}
		Ze2 += w3[p] * (0.25f * (rou_f ^ rou_sc) + 1.0f / (k_ * k_))*Ztemp;
	}

	Ze = (Ze1 + 1.0f / S_s * Ze2);// +Ze3;

	return Ze;
}

Complex SED_CM_FMM::MMZmpSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	/*Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0);// Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	VectorR3 vsc = (vs3[0] + vs3[1] + vs3[2]) / 3.0;
	VectorR3 rou_sc = vsc - vs; //.1.

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
	//电场  非奇异部分
	rou_f = vf - vgf[p];
	for (int q = 0; q < 3; q++)
	{
	rou_s = vgs[q] - vs; //.2.
	value_t r = (vgf[p] - vgs[q]).Norm();
	Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
	Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * coff; //.3.
	}
	//电场  奇异部分

	Complex Ztemp(0.0, 0.0);
	value_t Is2(0.0f);
	VectorR3 Is1(0.0f, 0.0f, 0.0f);
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[1]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vgf[p] - ni*(ni^vgf[p]);
	VectorR3 rou_n = vs - ni*(ni^vs);
	for (int i = 0; i < 3; i++)
	{
	VectorR3 li = vs3[(i + 1) % 3] - vs3[i];//not normalize
	VectorR3 lj = vs3[(i + 2) % 3] - vs3[(i + 1) % 3];//not normalize
	VectorR3 r_f = vgf[p] - vs3[i];
	VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
	VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
	VectorR3 u = li*n;
	VectorR3 ui = u / u.Norm();//normalize
	VectorR3 li_1 = li / li.Norm();//normalize
	value_t p0 = abs((vs3[i] - vgf[p]) ^ ui);
	value_t d = (vgf[p] - vs3[i]) ^ ni;
	value_t lp = (rou_p - rou) ^ li_1;
	value_t lm = (rou_m - rou) ^ li_1;
	value_t rp = sqrt(p0*p0 + d*d + lp*lp);
	value_t rm = sqrt(p0*p0 + d*d + lm*lm);

	if (p0 < 1e-10)
	continue;
	VectorR3 pi = ((rou_p - rou) - lp * li_1).Normalize();
	value_t r02 = p0*p0 + d*d;
	Is1 += (0.5f*ui)*(r02*log((rp + lp) / (rm + lm)) + (lp*rp) - (lm*rm));
	Is2 += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
	}
	Ze2 += w3[p] * (0.25f * (rou_f ^ (Is1 + Is2*(rou - rou_n))) + (Is2 / (k_ * k_)));
	}*/


	/*if (abs(Sigma) > 1e-8)
	{
	for (int i = 0; i < 3; i++)
	{
	rou_f = vgf[i] - vf;
	rou_s = vgf[i] - vs;
	Ze3 += w3[i] * rou_f^rou_s;
	}
	Ze3 = (PI / (S_s*J0*k_*Z0*Sigma))*Ze3;
	}*/
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	VectorR3 vsc = (vs3[0] + vs3[1] + vs3[2]) / 3.0;
	VectorR3 rou_sc = vsc - vs; //.1.

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
		//电场  非奇异部分
		rou_f = vf - vgf[p];
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs; //.2.
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * coff; //.3.
		}
		//电场  奇异部分

		Complex Ztemp(0.0, 0.0);
		for (int i = 0; i < 3; i++)
		{
			VectorR3 Li = vs3[(i + 1) % 3] - vs3[i];
			Li.Normalize();//边向量

			value_t Lpi = (vs3[(i + 1) % 3] - vgf[p]) ^ Li; //(边正顶点 - 场点)·(边向量)
			value_t Lmi = (vs3[i] - vgf[p]) ^ Li; //(边负顶点 - 场点)·(边向量)
			value_t Rpi = (vs3[(i + 1) % 3] - vgf[p]).Norm(); //边正顶点 到 场点 的距离
			value_t Rmi = (vs3[i % 3] - vgf[p]).Norm(); //边负顶点 到 场点 的距离
			VectorR3 sj0 = vs3[i] + Li * ((vgf[p] - vs3[i]) ^ Li); //垂足在u方向位置矢量
			value_t P0i = (sj0 - vgf[p]).Norm(); //垂足 到 场点 的距离
			if (P0i > 1e-10)
				Ztemp += P0i * log((Rpi + Lpi) / (Rmi + Lmi)); //.4.
		}
		Ze2 += w3[p] * (0.25f * (rou_f ^ rou_sc) + 1.0f / (k_ * k_))*Ztemp;
	}

	Ze = (Ze1 + 1.0f / S_s * Ze2);// +Ze3;

	return Ze;
}

Complex SED_CM_FMM::MMZmmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	/*Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0);// Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	VectorR3 vsc = (vs3[0] + vs3[1] + vs3[2]) / 3.0;
	VectorR3 rou_sc = vs - vsc; //.1.

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
	//电场  非奇异部分
	rou_f = vf - vgf[p];
	for (int q = 0; q < 3; q++)
	{
	rou_s = vs - vgs[q]; //.2.
	value_t r = (vgf[p] - vgs[q]).Norm();
	Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
	Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * coff; //.3.
	}
	//电场  奇异部分

	Complex Ztemp(0.0, 0.0);
	value_t Is2(0.0f);
	VectorR3 Is1(0.0f, 0.0f, 0.0f);
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[1]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vgf[p] - (ni^vgf[p])*ni;
	VectorR3 rou_n = vs - (ni^vs)*ni;
	for (int i = 0; i < 3; i++)
	{
	VectorR3 li = vs3[(i + 1) % 3] - vs3[i];//not normalize
	VectorR3 lj = vs3[(i + 2) % 3] - vs3[(i + 1) % 3];//not normalize
	VectorR3 r_f = vgf[p] - vs3[i];
	VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
	VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
	VectorR3 u = li*n;
	VectorR3 ui = u / u.Norm();//normalize
	VectorR3 li_1 = li / li.Norm();//normalize
	value_t p0 = abs((vs3[i] - vgf[p]) ^ ui);
	value_t d = (vgf[p] - vs3[i]) ^ ni;
	value_t lp = (rou_p - rou) ^ li_1;
	value_t lm = (rou_m - rou) ^ li_1;
	value_t rp = sqrt(p0*p0 + d*d + lp*lp);
	value_t rm = sqrt(p0*p0 + d*d + lm*lm);

	if (p0 < 1e-10)
	continue;
	VectorR3 pi = ((rou_p - rou) - lp * li_1).Normalize();
	value_t r02 = p0*p0 + d*d;
	Is1 += (0.5f*ui)*(r02*log((rp + lp) / (rm + lm)) + (lp*rp) - (lm*rm));
	Is2 += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
	}
	Ze2 += w3[p] * (0.25f * (-0.1f*rou_f ^ (Is1 + Is2*(rou - rou_n))) - (Is2 / (k_ * k_)));
	}*/

	/*if (abs(Sigma) > 1e-8)
	{
	for (int i = 0; i < 3; i++)
	{
	rou_f = vgf[i] - vf;
	rou_s = vgf[i] - vs;
	Ze3 += w3[i] * rou_f^rou_s;
	}
	Ze3 = (PI / (S_s*J0*k_*Z0*Sigma))*Ze3;
	}*/
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	VectorR3 vsc = (vs3[0] + vs3[1] + vs3[2]) / 3.0;
	VectorR3 rou_sc = vs - vsc; //.1.

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
		//电场  非奇异部分
		rou_f = vf - vgf[p];
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vgs[q]; //.2.
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * coff; //.3.
		}
		//电场  奇异部分

		Complex Ztemp(0.0, 0.0);
		for (int i = 0; i < 3; i++)
		{
			VectorR3 Li = vs3[(i + 1) % 3] - vs3[i];
			Li.Normalize();//边向量

			value_t Lpi = (vs3[(i + 1) % 3] - vgf[p]) ^ Li; //(边正顶点 - 场点)·(边向量)
			value_t Lmi = (vs3[i] - vgf[p]) ^ Li; //(边负顶点 - 场点)·(边向量)
			value_t Rpi = (vs3[(i + 1) % 3] - vgf[p]).Norm(); //边正顶点 到 场点 的距离
			value_t Rmi = (vs3[i % 3] - vgf[p]).Norm(); //边负顶点 到 场点 的距离
			VectorR3 sj0 = vs3[i] + Li * ((vgf[p] - vs3[i]) ^ Li); //垂足在u方向位置矢量
			value_t P0i = (sj0 - vgf[p]).Norm(); //垂足 到 场点 的距离
			if (P0i > 1e-10)
				Ztemp += P0i * log((Rpi + Lpi) / (Rmi + Lmi)); //.4.
		}
		Ze2 += w3[p] * (0.25f * (rou_f ^ rou_sc) - 1.0f / (k_ * k_))*Ztemp;
	}

	Ze = (Ze1 + 1.0f / S_s * Ze2);// +Ze3;

	return Ze;
}

Complex SED_CM_FMM::MMEDMKernel(VectorR3 *vf_p3, VectorR3 *vf_m3, VectorR3 *vs_p3, VectorR3 *vs_m3, value_t &l_fld, value_t &l_src)
{
	VectorR3 vf_cen_plu = centerv3(vf_p3);
	VectorR3 vf_cen_min = centerv3(vf_m3);

	VectorR3 vs_cen_plu = centerv3(vs_p3);
	VectorR3 vs_cen_min = centerv3(vs_m3);

	VectorR3 r_fld = (vf_cen_plu + vf_cen_min) / 2.0f;
	VectorR3 r_src = (vs_cen_plu + vs_cen_min) / 2.0f;

	VectorR3 R = r_fld - r_src;
	value_t R_norm = R.Norm();

	VectorR3 m_s = l_src * (vs_cen_min - vs_cen_plu);
	Complex C = (1.0f / (R_norm*R_norm))*(1.0f + 1.0f / (J0*k_*R_norm));
	VectorR3 M_s = (R^m_s)*R / (R_norm*R_norm);

	VectorC3 E;
	E = ((M_s - m_s)*((J0*k_ / R_norm) + C) + 2.0f*M_s*C)*exp(-J0 * k_*R_norm)*Z0 / PI4;

	return (-l_fld * (E ^ (vf_cen_min - vf_cen_plu)));
}

value_t SED_CM_FMM::IsSingular(VectorR3 * vs3, VectorR3 & vf)
{
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[1]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vf - ni * (ni^vf);

	value_t Is = 0.0f;
	value_t d = (vf - vs3[0]) ^ ni;
	for (int i = 0; i < 3; i++)
	{
		VectorR3 li = (vs3[(i + 1) % 3] - vs3[i]).Normalize();//not normalize
		VectorR3 r_f = vf - vs3[i];
		VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
		VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
		VectorR3 ui = li * ni;
		value_t p0 = abs((vs3[i] - vf) ^ ui);
		value_t lp = (vs3[(i + 1) % 3] - vf) ^ li;
		value_t lm = (vs3[i] - vf) ^ li;
		value_t lp1 = (vf - vs3[(i + 1) % 3]) ^ li;
		value_t lm1 = (vf - vs3[i]) ^ li;
		value_t rp = sqrt(p0*p0 + d * d + lp * lp);
		value_t rm = sqrt(p0*p0 + d * d + lm * lm);

		if (p0 < 1e-8)
			continue;

		VectorR3 pi = ((rou_p - rou) - lp * li).Normalize();
		value_t r02 = p0 * p0 + d * d;
		Is += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
	}
	return Is;
}

VectorR3 SED_CM_FMM::IspSingular(VectorR3 * vs3, VectorR3 & vf)
{
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[1]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vf - ni * (ni^vf);

	VectorR3 Isp(0, 0, 0);
	value_t d = (vf - vs3[0]) ^ ni;
	for (int i = 0; i < 3; i++)
	{
		VectorR3 li = (vs3[(i + 1) % 3] - vs3[i]).Normalize();//not normalize
		VectorR3 r_f = vf - vs3[i];
		VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
		VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
		VectorR3 ui = li * ni;
		value_t p0 = abs((vs3[i] - vf) ^ ui);
		value_t lp = (vs3[(i + 1) % 3] - vf) ^ li;
		value_t lm = (vs3[i] - vf) ^ li;
		value_t lp1 = (vf - vs3[(i + 1) % 3]) ^ li;
		value_t lm1 = (vf - vs3[i]) ^ li;
		value_t rp = sqrt(p0*p0 + d * d + lp * lp);
		value_t rm = sqrt(p0*p0 + d * d + lm * lm);

		if (p0 < 1e-8)
			continue;

		VectorR3 pi = ((rou_p - rou) - lp * li).Normalize();
		value_t r02 = p0 * p0 + d * d;
		//Isp += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
		Isp += (r02*(log((rp + lp) / (rm + lm))) + rp * lp - rm * lm)*ui;
	}
	return 0.5f*Isp;
}

value_t SED_CM_FMM::IvSingular(VectorR3 * vs3, VectorR3 & vf)
{
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[0]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vf - ni * (ni^vf);

	value_t Iv = 0.0f;
	value_t d = (vf - vs3[0]) ^ ni;
	for (int i = 0; i < 3; i++)
	{
		VectorR3 li = (vs3[(i + 1) % 3] - vs3[i]).Normalize();//not normalize
		VectorR3 r_f = vf - vs3[i];
		VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
		VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
		VectorR3 ui = li * ni;
		value_t p0 = abs((vs3[i] - vf) ^ ui);
		value_t lp = (vs3[(i + 1) % 3] - vf) ^ li;
		value_t lm = (vs3[i] - vf) ^ li;
		value_t rp = sqrt(p0*p0 + d * d + lp * lp);
		value_t rm = sqrt(p0*p0 + d * d + lm * lm);

		if (p0 < 1e-8)
			continue;

		VectorR3 pi = ((rou_p - rou) - lp * li).Normalize();
		value_t r02 = p0 * p0 + d * d;
		Iv += (pi^ui)*(abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))) - p0 * log((rp + lp) / (rm + lm)));
	}
	return d * Iv;
}

VectorR3 SED_CM_FMM::IvpSingular(VectorR3 * vs3, VectorR3 & vf)
{
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[0]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vf - ni * (ni^vf);
	value_t Ivp(0.0f);
	value_t d = (vf - vs3[0]) ^ ni;
	for (int i = 0; i < 3; i++)
	{
		VectorR3 li = (vs3[(i + 1) % 3] - vs3[i]).Normalize();//not normalize

		VectorR3 r_f = vf - vs3[i];
		VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
		VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
		VectorR3 ui = li * ni;
		value_t p0 = abs((vs3[i] - vf) ^ ui);
		value_t lp = (vs3[(i + 1) % 3] - vf) ^ li;
		value_t lm = (vs3[i] - vf) ^ li;
		value_t rp = sqrt(p0*p0 + d * d + lp * lp);
		value_t rm = sqrt(p0*p0 + d * d + lm * lm);

		if (p0 < 1e-8)
			continue;

		VectorR3 pi = ((rou_p - rou) - lp * li).Normalize();
		value_t r02 = p0 * p0 + d * d;
		value_t d3 = abs(d*d*d);
		value_t d2 = d * d;
		Ivp += (pi^ui)*(((p0*(r02 + 2.0f*d2)) / 2.0f)*log((rp + lp) / (rm + lm)) + (p0 / 2.0f)*(lp*rp - lm * rm) - d3 * (atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
	}
	return Ivp * ni;
}

void SED_CM_FMM::TetToTri(VectorR3 * v4, VectorR3 ** v4_3)
{
	v4_3[0][0] = v4[0];
	v4_3[0][1] = v4[2];
	v4_3[0][2] = v4[1];

	v4_3[1][0] = v4[0];
	v4_3[1][1] = v4[3];
	v4_3[1][2] = v4[2];

	v4_3[2][0] = v4[0];
	v4_3[2][1] = v4[1];
	v4_3[2][2] = v4[3];

	v4_3[3][0] = v4[1];
	v4_3[3][1] = v4[2];
	v4_3[3][2] = v4[3];
}

Complex SED_CM_FMM::z1(VectorR3 * vf4, VectorR3 * vs4, VectorR3 & vf, VectorR3 & vs, value_t & epsp, value_t & epsm)
{

	return Complex();
}

Complex SED_CM_FMM::DVKernel(VectorR3 *vp4, VectorR3 *vm4, VectorR3 &vp, VectorR3 &vm, int &_nvm)
{
	VectorR3 vgp[4], vgm[4];
	VectorR3 rou_p, rou_m; //ρ+-
	Complex Ve;
	value_t kp = 0.0f, km = 0.0f;
	//kp = (epsp - 1.0f) / epsp;
	//km = (epsm - 1.0f) / epsm;

	Gauss4Point(vp4[0], vp4[1], vp4[2], vp4[3], vgp);
	Gauss4Point(vm4[0], vm4[1], vm4[2], vm4[3], vgm);

	VectorC3 Ei(0, 0, 0);
	Complex Vep(0, 0), Vem(0, 0), G(0, 0);

	for (int a = 0; a < 4; a++)
	{
		//////////////////////////////////////////////////////////////////////////
		rou_p = vgp[a] - vp;
		G = exp(-J0 * k_ * (vgp[a] ^ inc_k_));

		Ei = G * inc_e_;
		Vep += wt4[a] * (rou_p ^ Ei);

		//////////////////////////////////////////////////////////////////////////
		if (_nvm != -1)
		{
			rou_m = vm - vgm[a];
			G = exp(-J0 * k_ * (vgm[a] ^ inc_k_));

			Ei = G * inc_e_;
			Vem += wt4[a] * (rou_m ^ Ei);
		}

	}

	Ve = Vep + Vem;

	return Ve;

}

Complex SED_CM_FMM::MVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &vp, VectorR3 &vm)
{
	VectorR3 vgp[7], vgm[7];
	VectorR3 rou_p, rou_m; //ρ+-
	Complex Ve;

	Gauss7Point(vp3[0], vp3[1], vp3[2], vgp);
	Gauss7Point(vm3[0], vm3[1], vm3[2], vgm);

	VectorC3 Ei(0, 0, 0)/*, Hi(0, 0, 0)*/;
	Complex Vep(0, 0), Vem(0, 0)/*, Vmp(0, 0), Vmm(0, 0)*/, G(0, 0);

	for (int a = 0; a < 7; a++)
	{
		//////////////////////////////////////////////////////////////////////////
		rou_p = vgp[a] - vp;
		G = exp(-J0 * k_ * (vgp[a] ^ inc_k_));

		Ei = G * inc_e_;
		Vep += w7[a] * (rou_p ^ Ei);

		//////////////////////////////////////////////////////////////////////////
		rou_m = vm - vgm[a];
		G = exp(-J0 * k_ * (vgm[a] ^ inc_k_));

		Ei = G * inc_e_;
		Vem += w7[a] * (rou_m ^ Ei);
	}

	Ve = (Vep + Vem);

	return Ve;
	//后面要乘以 （0.5*L）
}

Complex SED_CM_FMM::FSPGF(VectorR3 &r, int _t_sum, value_t _Dx, value_t _Dy)
{
	Complex PG(0.0, 0.0), PG1(0.0, 0.0), PG2(0.0, 0.0);
	VectorR3 R = r;
	value_t H;
	double delta = 0.0;
	Complex a(0.0, 0.0);
	H = sqrt(PI / (Dx*Dy));
	for (int i = -_t_sum; i <= _t_sum; i++)
	{
		for (int j = -_t_sum; j <= _t_sum; j++)
		{

			delta = ((float)i*PI / Dx)*((float)i*PI / Dx) + ((float)j*PI / Dy) *((float)j*PI / Dy);
			//R = ((d_x - m*Dx) ^ 2 + (d_y - n*Dy) ^ 2 + (d_z) ^ 2);
			R.x = r.x - i * Dx;
			R.y = r.y - j * Dy;
			value_t R_N = R.Norm();

			if (delta >= ((k_ * k_) / 4.0f))
				a = sqrt(delta - ((k_ * k_) / 4.0f));
			else
				a = J0 * (value_t)sqrt(((k_ * k_) / 4.0f - delta));

			PG1 += ((exp(-2.0f*J0 * PI*(i*r.x / Dx + j * r.y / Dy))) / a)*((exp(2.0f * a*r.z))*(erfc(r.z*H + a / H)) + (exp(-2.0f * a*r.z))*(erfc((a / H) - r.z*H)));
			PG2 += (exp(-J0 * k_*(inc_k_.x*i*Dx + inc_k_.y*j*Dy)) / R_N)*(erfc(R_N*H + (-J0 * k_) / (2.0f * H))*exp(-J0 * k_*R_N)).real();
			/*if (isinf(PG1.real()) || isinf(PG1.imag()) || isinf(PG2.imag()) || isinf(PG2.real()))
			{
			Qcout << "This number isn't a number!" << std::endl;
			Qcout << "please cheak this part!" << std::endl;
			}*/
			if (_isnan(PG1.real()) || _isnan(PG2.real()))
			{
				Qcout << "This element is not a number!" << std::endl;
				Qcout << "Please check this part!" << std::endl;
			}
		}
	}

	PG = (PI*PG1) / (2.0f * Dx*Dy) + PG2;
	return PG;
}

Complex SED_CM_FMM::erfz(Complex z)
{
	/*std::complex<double> j0(0.0, 1.0);
	double pi = 3.1415926;
	value_t sqrtpi = sqrt(pi);
	value_t abs_z = abs(z);



	Complex f;

	if (abs_z <= 8)
	{
	int n = 32;
	double x = z.real(), y = z.imag();
	std::complex<double> k1 = 2 / pi * exp(-x*x), k2 = exp(-2.0 * j0  * x*y);
	std::complex<double> s1 = erf(x);
	std::complex<double> s2(0.0, 0.0), s3(0.0, 0.0), s4(0.0, 0.0), s5(0.0, 0.0), s6(0.0, 0.0);
	if (abs(x) > 1e-5)
	{
	s2 = k1 / (4.0 * x)*(1.0 - k2);
	}
	else
	{
	s2 = (j0 / pi) *y;
	}

	f = s1 + s2;

	double xk=0.0, yk=0.0;

	if (abs(y) > 1e-5)
	{
	xk = x;
	yk = y;
	}

	s5 = 0;
	for (int i = 1; i <= n; i++)
	{
	s3 = exp(-(1.0*(i*i)) / 4.0) / (1.0*(i*i + 4 * xk*xk));
	s4 = 2.0 * xk - k2*(2.0 * xk*cosh(i*yk) - 1.0*i*j0*sinh(i*yk));
	s5 = s5 + s3*s4;
	}
	s6 = k1*s5;
	f = f + (arma::cx_float)s6;

	return f;
	}
	else
	{
	if (z.real() < 0)
	z = -z;
	int nmax = 193;
	Complex s = (1.0, 0.0);
	Complex y = 2.0f * z*z;
	for (int i = nmax; i >= 1; i -= 2)
	{
	s = 1.0f - 1.0f*i*(s / y);
	}

	f = 1.0f - s*exp(-z*z) / (sqrtpi*z);

	if (z.real() < 0)
	f = -f;

	if (abs(z.real()) < 1e-8)
	f = f - 1.0f;

	return f;
	}*/

	double x = z.real(), y = z.imag();
	std::complex<double> k1 = 2.0 / PI * exp(-x * x), k2 = exp(-2.0 * (std::complex<double>)J0  * x*y);
	std::complex<double> s1 = erf(x);
	std::complex<double> s2(0.0, 0.0), s3(0.0, 0.0), s4(0.0, 0.0), s5(0.0, 0.0), s6(0.0, 0.0);
	int n = 25;
	Complex f;
	if (abs(x) > 1e-5)
	{
		s2 = k1 / (4.0 * x)*(1.0 - k2);
	}
	else
	{
		s2 = ((std::complex<double>)J0 / (double)PI) *y;
	}


	s5 = 0;
	for (int i = 1; i <= n; i++)
	{
		s3 = exp(-(1.0*(i*i)) / 4.0) / (1.0*(i*i + 4.0 * x*x));
		s4 = 2.0 * x - k2 * (2.0 * x*cosh(i*y) - 1.0*i*(std::complex<double>)J0*sinh(i*y));
		s5 = s5 + s3 * s4;
	}
	s6 = k1 * s5;
	f = (arma::cx_float)(s1 + s2 + s6);

	return f;
}

void SED_CM_FMM::KernelIntegral()
{
	Tetrahedron tet_i, tet_j, tet_m, tet_n;
	VectorR3 v_tet_i[4], v_tet_j[4], v_tet_m[4], v_tet_n[4];
	VectorR3 vgi[4], vgj[4], vgm[4], vgn[4];

	Qcout << '\n' << "Fill the kernelintegral:";
	tool::BarAndPercent bar;
	for (int i = 0; i < 9; i++)
	{
		for (int j = i + 1; j < 9; j++)
		{
			for (int m = 0; m < tet_num; m++)
			{
				tet_m = mesh_ptr_->getTetrahedronRef(m);
				tet_m.GetVertices(v_tet_m[0], v_tet_m[1], v_tet_m[2], v_tet_m[3]);
				shiftcoordv4(i, v_tet_m);
				Gauss4Point(v_tet_m[0], v_tet_m[1], v_tet_m[2], v_tet_m[3], vgm);
				for (int n = 0; n < tet_num; n++)
				{
					tet_n = mesh_ptr_->getTetrahedronRef(n);
					tet_n.GetVertices(v_tet_n[0], v_tet_n[1], v_tet_n[2], v_tet_n[3]);
					shiftcoordv4(j, v_tet_n);
					Gauss4Point(v_tet_n[0], v_tet_n[1], v_tet_n[2], v_tet_n[3], vgn);
					FillZKernel(vgm, vgn, m, n, v_tet_n);
				}
			}
		}
	}
	for (int i = 0; i < tet_num; i++)
	{
		bar(i + 1, tet_num);
		tet_i = mesh_ptr_->getTetrahedronRef(i);
		tet_i.GetVertices(v_tet_i[0], v_tet_i[1], v_tet_i[2], v_tet_i[3]);
		Gauss4Point(v_tet_i[0], v_tet_i[1], v_tet_i[2], v_tet_i[3], vgi);
		for (int j = i; j < tet_num; j++)
		{
			tet_j = mesh_ptr_->getTetrahedronRef(j);
			tet_j.GetVertices(v_tet_j[0], v_tet_j[1], v_tet_j[2], v_tet_j[3]);
			Gauss4Point(v_tet_j[0], v_tet_j[1], v_tet_j[2], v_tet_j[3], vgj);
			FillZKernel(vgi, vgj, i, j, v_tet_j);
		}
	}
	Qcout << "test" << std::endl;
}

SED_CM_FMM::KIM SED_CM_FMM::SED_KernelIntegral(int &sedid_f, int &sedid_s)
{
	Tetrahedron tet_i, tet_j, tet_m, tet_n;
	VectorR3 v_tet_i[4], v_tet_j[4], v_tet_m[4], v_tet_n[4];
	VectorR3 vgi[4], vgj[4], vgm[4], vgn[4];
	KIM ik(tet_num);

	if (sedid_f == sedid_s)
	{
		for (int i = 0; i < tet_num; i++)
		{
			tet_i = mesh_ptr_->getTetrahedronRef(i);
			tet_i.GetVertices(v_tet_i[0], v_tet_i[1], v_tet_i[2], v_tet_i[3]);
			shiftcoordv4(sedid_f, v_tet_i);
			Gauss4Point(v_tet_i[0], v_tet_i[1], v_tet_i[2], v_tet_i[3], vgi);
			for (int j = i; j < tet_num; j++)
			{
				tet_j = mesh_ptr_->getTetrahedronRef(j);
				tet_j.GetVertices(v_tet_j[0], v_tet_j[1], v_tet_j[2], v_tet_j[3]);
				shiftcoordv4(sedid_s, v_tet_j);
				Gauss4Point(v_tet_j[0], v_tet_j[1], v_tet_j[2], v_tet_j[3], vgj);
				FillZKernel_SED_COM(vgi, vgj, i, j, v_tet_j, ik);
			}
		}
	}
	else
	{
		for (int m = 0; m < tet_num; m++)
		{
			tet_m = mesh_ptr_->getTetrahedronRef(m);
			tet_m.GetVertices(v_tet_m[0], v_tet_m[1], v_tet_m[2], v_tet_m[3]);
			shiftcoordv4(sedid_f, v_tet_m);
			Gauss4Point(v_tet_m[0], v_tet_m[1], v_tet_m[2], v_tet_m[3], vgm);
			for (int n = 0; n < tet_num; n++)
			{
				tet_n = mesh_ptr_->getTetrahedronRef(n);
				tet_n.GetVertices(v_tet_n[0], v_tet_n[1], v_tet_n[2], v_tet_n[3]);
				shiftcoordv4(sedid_s, v_tet_n);
				Gauss4Point(v_tet_n[0], v_tet_n[1], v_tet_n[2], v_tet_n[3], vgn);
				FillZKernel_SED_DIFF(vgm, vgn, m, n, v_tet_n, ik);
			}
		}
	}

	return ik;
}

SED_CM_FMM::KIM SED_CM_FMM::RED_KernelIntegral(int &xid_f, int &yid_f, int &xid_s, int &yid_s)
{
	Tetrahedron tet_i, tet_j, tet_m, tet_n;
	VectorR3 v_tet_i[4], v_tet_j[4], v_tet_m[4], v_tet_n[4];
	VectorR3 vgi[4], vgj[4], vgm[4], vgn[4];
	KIM ik(tet_num);
	//int arrid_f = ArrayToSED(xid_f, yid_f);
	//int arrid_s = ArrayToSED(xid_s, yid_s);
	if (xid_f == xid_s && yid_f == yid_s)
	{
		for (int i = 0; i < tet_num; i++)
		{
			tet_i = mesh_ptr_->getTetrahedronRef(i);
			tet_i.GetVertices(v_tet_i[0], v_tet_i[1], v_tet_i[2], v_tet_i[3]);
			shiftcoordv4(xid_f, yid_f, v_tet_i);
			Gauss4Point(v_tet_i[0], v_tet_i[1], v_tet_i[2], v_tet_i[3], vgi);
			for (int j = i; j < tet_num; j++)
			{
				tet_j = mesh_ptr_->getTetrahedronRef(j);
				tet_j.GetVertices(v_tet_j[0], v_tet_j[1], v_tet_j[2], v_tet_j[3]);
				shiftcoordv4(xid_s, yid_s, v_tet_j);
				Gauss4Point(v_tet_j[0], v_tet_j[1], v_tet_j[2], v_tet_j[3], vgj);

				FillZKernel_SED_COM(vgi, vgj, i, j, v_tet_j, ik);
			}
		}
	}
	else
	{
		for (int m = 0; m < tet_num; m++)
		{
			tet_m = mesh_ptr_->getTetrahedronRef(m);
			tet_m.GetVertices(v_tet_m[0], v_tet_m[1], v_tet_m[2], v_tet_m[3]);
			shiftcoordv4(xid_f, yid_f, v_tet_m);
			Gauss4Point(v_tet_m[0], v_tet_m[1], v_tet_m[2], v_tet_m[3], vgm);
			for (int n = 0; n < tet_num; n++)
			{
				tet_n = mesh_ptr_->getTetrahedronRef(n);
				tet_n.GetVertices(v_tet_n[0], v_tet_n[1], v_tet_n[2], v_tet_n[3]);
				shiftcoordv4(xid_s, yid_s, v_tet_n);
				Gauss4Point(v_tet_n[0], v_tet_n[1], v_tet_n[2], v_tet_n[3], vgn);
				FillZKernel_SED_DIFF(vgm, vgn, m, n, v_tet_n, ik);
			}
		}
	}

	return ik;
}

void SED_CM_FMM::FillZKernel(VectorR3 *vm4, VectorR3 *vn4, int &m, int &n, VectorR3 *v4)
{
	Complex G(0, 0);
	value_t r;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			r = (vm4[i] - vn4[j]).Norm();

			if (m != n)
			{
				G = exp(-J0 * k_*r) / r;

				Z2K1(m, n) += wt4[i] * wt4[j] * (vm4[i] ^ vn4[j])*G;
				Z2K1(n, m) = Z2K1(m, n);

				Z2K2(m, n) += wt4[i] * wt4[j] * vn4[j] * G;
				Z2K3(m, n) += wt4[i] * wt4[j] * vm4[i] * G;

				Z2K2(n, m) = Z2K3(m, n);
				Z2K3(n, m) = Z2K2(m, n);

				Z2K4(m, n) += wt4[i] * wt4[j] * G;
				Z2K4(n, m) = Z2K4(m, n);

				Z5K(m, n) = Z2K4(m, n);
				Z5K(n, m) = Z2K4(m, n);
			}
			else
			{
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);

				Z2K1(m, m) += wt4[i] * wt4[j] * (vm4[i] ^ vn4[j])*G;
				Z2K2(m, m) += wt4[i] * wt4[j] * vn4[j] * G;
				Z2K3(m, m) += wt4[i] * wt4[j] * vm4[i] * G;
				Z2K4(m, m) += wt4[i] * wt4[j] * G;
				Z5K(m, m) = Z2K4(m, m);
			}

		}
		if (m == n)
		{
			VectorR3 Ivp(0.0f, 0.0f, 0.0f);
			value_t Iv = 0.0f;

			Pv4_3[0] = v4_3[0];
			Pv4_3[1] = v4_3[1];
			Pv4_3[2] = v4_3[2];
			Pv4_3[3] = v4_3[3];
			TetToTri(v4, Pv4_3);

			for (int k = 0; k < 4; k++)
			{
				Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[k], vm4[i]);
				Iv += 0.5f*IvSingular(v4_3[k], vm4[i]);
			}
			RIVP(m) += wt4[i] * (vm4[i] ^ Ivp);
			RRIV(m) += wt4[i] * (vm4[i] ^ vm4[i])*Iv;
			RIV[m] += wt4[i] * Iv*vm4[i];
			IVP[m] += wt4[i] * Ivp;
			IV(m) += wt4[i] * Iv;
		}
	}
}

void SED_CM_FMM::FillZKernel_SED_COM(VectorR3 *vm4, VectorR3 *vn4, int &m, int &n, VectorR3 *v4, KIM& ik)
{
	Complex G(0, 0);
	value_t r;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			r = (vm4[i] - vn4[j]).Norm();

			if (m != n)
			{
				G = exp(-J0 * k_*r) / r;

				ik.Z2K1(m, n) += wt4[i] * wt4[j] * (vm4[i] ^ vn4[j])*G;
				ik.Z2K1(n, m) = ik.Z2K1(m, n);

				ik.Z2K2(m, n) += wt4[i] * wt4[j] * vn4[j] * G;
				ik.Z2K3(m, n) += wt4[i] * wt4[j] * vm4[i] * G;

				ik.Z2K2(n, m) = ik.Z2K3(m, n);
				ik.Z2K3(n, m) = ik.Z2K2(m, n);

				ik.Z2K4(m, n) += wt4[i] * wt4[j] * G;
				ik.Z2K4(n, m) = ik.Z2K4(m, n);

				ik.Z5K(m, n) = ik.Z2K4(m, n);
				ik.Z5K(n, m) = ik.Z2K4(m, n);
			}
			else
			{
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);

				ik.Z2K1(m, m) += wt4[i] * wt4[j] * (vm4[i] ^ vn4[j])*G;
				ik.Z2K2(m, m) += wt4[i] * wt4[j] * vn4[j] * G;
				ik.Z2K3(m, m) += wt4[i] * wt4[j] * vm4[i] * G;
				ik.Z2K4(m, m) += wt4[i] * wt4[j] * G;
				ik.Z5K(m, m) = ik.Z2K4(m, m);
			}

		}
		if (m == n)
		{
			VectorR3 Ivp(0.0f, 0.0f, 0.0f);
			value_t Iv = 0.0f;

			Pv4_3[0] = v4_3[0];
			Pv4_3[1] = v4_3[1];
			Pv4_3[2] = v4_3[2];
			Pv4_3[3] = v4_3[3];
			TetToTri(v4, Pv4_3);

			for (int k = 0; k < 4; k++)
			{
				Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[k], vm4[i]);
				Iv += 0.5f*IvSingular(v4_3[k], vm4[i]);
			}
			ik.RIVP(m) += wt4[i] * (vm4[i] ^ Ivp);
			ik.RRIV(m) += wt4[i] * (vm4[i] ^ vm4[i])*Iv;
			ik.RIV[m] += wt4[i] * Iv*vm4[i];
			ik.IVP[m] += wt4[i] * Ivp;
			ik.IV(m) += wt4[i] * Iv;
		}
	}
}

void SED_CM_FMM::FillZKernel_SED_DIFF(VectorR3 *vm4, VectorR3 *vn4, int &m, int &n, VectorR3 *v4, KIM& ik)
{
	Complex G(0, 0);
	value_t r;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			r = (vm4[i] - vn4[j]).Norm();
			G = exp(-J0 * k_*r) / r;

			ik.Z2K1(m, n) += wt4[i] * wt4[j] * (vm4[i] ^ vn4[j])*G;
			ik.Z2K2(m, n) += wt4[i] * wt4[j] * vn4[j] * G;
			ik.Z2K3(m, n) += wt4[i] * wt4[j] * vm4[i] * G;
			ik.Z2K4(m, n) += wt4[i] * wt4[j] * G;

			ik.Z5K(m, n) = ik.Z2K4(m, n);
		}

	}
}


void SED_CM_FMM::CalculatePGFgrid()
{
	grid_x = std::ceil((max_boundary.x - min_boundary.x) / d_x) + 2;
	grid_y = std::ceil((max_boundary.y - min_boundary.y) / d_y) + 2;
	grid_z = std::ceil((max_boundary.z - min_boundary.z) / d_z) + 2;

	//int grid_t = grid_x*grid_y*grid_z;
	tool::BarAndPercent bar;
	for (int z = 0; z < grid_z; z++)
	{
		bar(z + 1, grid_z);
		for (int y = 0; y < grid_y; y++)
		{
			for (int x = 0; x < grid_x; x++)
			{
				if ((x == 0) && (y == 0) && (z == 0))
				{
					PGFgrid.push_back((0.0f, 0.0f));
					continue;
				}
				VectorR3 r;
				r.x = x * d_x;
				r.y = y * d_y;
				r.z = z * d_z;
				PGFgrid.push_back(FSPGF(r, t_sum, Dx, Dy));
			}
		}
	}
	Qcout << "The PGFgrid's filling is accomplished!" << std::endl;
}

Complex SED_CM_FMM::InterpolarPGF(VectorR3 &r)
{
	int x_seq, y_seq, z_seq;
	x_seq = floor(abs(r.x) / d_x);
	y_seq = floor(abs(r.y) / d_y);
	z_seq = floor(abs(r.z) / d_z);

	int n000 = z_seq * grid_x*grid_y + y_seq * grid_x + x_seq;
	int n100 = n000 + 1;
	int n010 = n000 + grid_x;
	int n110 = n010 + 1;
	int n001 = n000 + grid_x * grid_y;
	int n101 = n001 + 1;
	int n011 = n001 + grid_x;
	int n111 = n011 + 1;

	Complex G000 = PGFgrid[n000];
	Complex G100 = PGFgrid[n100];
	Complex G010 = PGFgrid[n010];
	Complex G110 = PGFgrid[n110];
	Complex G001 = PGFgrid[n001];
	Complex G101 = PGFgrid[n101];
	Complex G011 = PGFgrid[n011];
	Complex G111 = PGFgrid[n111];

	value_t effx = (abs(r.x) - x_seq * d_x) / d_x;
	value_t effy = (abs(r.y) - y_seq * d_y) / d_y;
	value_t effz = (abs(r.z) - z_seq * d_z) / d_z;

	Complex G = (1.0f - effx)*(1.0f - effy)*(1.0f - effz)*G000
		+ (1.0f - effx)*(1.0f - effy)*effz*G001
		+ (1.0f - effx)*effy*(1.0f - effz)*G010
		+ (1.0f - effx)*effy*effz*G011
		+ effx * (1.0f - effy)*(1.0f - effz)*G100
		+ effx * (1.0f - effy)*effz*G101
		+ effx * effy*(1.0f - effz)*G110
		+ effx * effy*effz*G111;

	return G;
}

void SED_CM_FMM::fillZ()
{
	KernelIntegral();
}

void SED_CM_FMM::fillZDD(int &start_f, int &end_f, int &start_s, int &end_s)
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3];
	Tetrahedron tet_plu, tet_min;
	value_t a_fld, a_src;//场与源的公共面面积
	value_t epsp, epsm;//源处介电常数
	int sedid_f, bedge_f, sedid_s, bedge_s;
	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp_f, vp_s, vm_f, vm_s, v_ct1, v_ct2, v_ct3;

	Complex coef(0.0f, 0.0f);
	Qcout << '\n' << "Fill the matrix ZDD:";
	//填充ZDD矩阵
	tool::BarAndPercent bar_perc1;
	coef = Z0 * Z0 / (J0* PI4);
	for (int f = start_f; f < end_f; f++)
	{
		bar_perc1(f + 1, unknowns_t);
		ct_ptr_->getBEDGEandSEDID(f, bedge_f, sedid_f);
		ct_ptr_->getCommonTriangle(f, vp_f, vm_f, t_fld_plu, t_fld_min, a_fld);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);
		if (vm_f != -1)
		{
			v_fld_plu = mesh_ptr_->getVertex(vp_f);//场顶点+
			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-

			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-

			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);
		}
		else
		{
			v_fld_plu = mesh_ptr_->getVertex(vp_f);//场点+
			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);
		}

		for (int s = start_s; s < end_s; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp_s, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);
			if (vm_s != -1)
			{
				v_src_plu = mesh_ptr_->getVertex(vp_s);
				v_src_min = mesh_ptr_->getVertex(vm_s);

				tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);

				tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);
			}
			else
			{
				v_src_plu = mesh_ptr_->getVertex(vp_s);
				tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
				tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);
			}

			//field tetrahedron+ <-->source tetrahedron+ 
			if (t_fld_plu == t_src_plu)
			{
				Zpp = DDZppSingular(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			}
			else
			{
				Zpp = DDZppKernel(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			}
			//field tetrahedron+ <-->source tetrahedron-
			if (t_fld_plu == t_src_min)
			{
				Zpm = DDZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			}
			else if (t_src_min != -1)
			{
				Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			}
			//field tetrahedron- <-->source tetrahedron+
			if (t_fld_min == t_src_plu)
			{
				Zmp = DDZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			}
			else if (t_fld_min != -1)
			{
				Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			}
			//field tetrahedron- <-->source tetrahedron-
			if ((t_fld_min != -1) && (t_src_min != -1))
			{
				if (t_fld_min == t_src_min)
				{
					Zmm = DDZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
				}
				else
					Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			}

			Z(f, s) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void SED_CM_FMM::fillZDM(int &start_f, int &end_f, int &start_s, int &end_s)
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle 
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	value_t a_fld, a_src;//SWG场与源的公共面面积
	value_t l_fld, l_src;//RWG场与源的公共边长度
	value_t epsp, epsm;//源处介电常数

					   //Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm_f, vm_s, vm, v_ct1, v_ct2, v_ct3;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f, 0.0f);
	Qcout << '\n' << "Fill the matrix ZMD:";
	//填充ZMD矩阵
	tool::BarAndPercent bar_perc2;
	coef = J0 * k_*Z0 / PI4;
	for (int f = 0; f < unknowns_t; f++)
	{
		bar_perc2(f + 1, unknowns_t);
		ct_ptr_->getCommonTriangle(f, vp, vm_f, t_fld_plu, t_fld_min, a_fld);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);
		if (vm_f != -1)
		{
			v_fld_plu = mesh_ptr_->getVertex(vp);//场顶点+
			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-

			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-

			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);
		}
		else
		{
			v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);
		}
		for (int s = 0; s < unknowns_e; s++)
		{
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);
			tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
			tri_min.getVertex(tri_min1, tri_min2, tri_min3);
			//fill
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);

			if (vm_f != -1)
			{
				Zmp = DMZmpKernel(v_fld_min4, v_src_plu3, v_fld_min, v_src_plu, v_fld_ct3, vm_f);
				Zmm = DMZmmKernel(v_fld_min4, v_src_min3, v_fld_min, v_src_min, v_fld_ct3, vm_f);
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else if ((tri_plu1 == v_ct1 && tri_plu2 == v_ct2 && tri_plu3 == v_ct3)
				|| (tri_plu1 == v_ct1 && tri_plu2 == v_ct3 && tri_plu3 == v_ct2)
				|| (tri_plu1 == v_ct2 && tri_plu2 == v_ct1 && tri_plu3 == v_ct3)
				|| (tri_plu1 == v_ct2 && tri_plu2 == v_ct3 && tri_plu3 == v_ct1)
				|| (tri_plu1 == v_ct3 && tri_plu2 == v_ct1 && tri_plu3 == v_ct2)
				|| (tri_plu1 == v_ct3 && tri_plu2 == v_ct2 && tri_plu3 == v_ct1))
			{
				Zpp = DMZppSingular(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
			}
			else if ((tri_min1 == v_ct1 && tri_min2 == v_ct2 && tri_min3 == v_ct3)
				|| (tri_min1 == v_ct1 && tri_min2 == v_ct3 && tri_min3 == v_ct2)
				|| (tri_min1 == v_ct2 && tri_min2 == v_ct1 && tri_min3 == v_ct3)
				|| (tri_min1 == v_ct2 && tri_min2 == v_ct3 && tri_min3 == v_ct1)
				|| (tri_min1 == v_ct3 && tri_min2 == v_ct1 && tri_min3 == v_ct2)
				|| (tri_min1 == v_ct3 && tri_min2 == v_ct2 && tri_min3 == v_ct1))
			{
				Zpm = DMZpmSingular(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else
			{
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}

			Z(f, s + unknowns_t) = coef * a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
		}

	}
}

void SED_CM_FMM::fillZMD(int &start_f, int &end_f, int &start_s, int &end_s)
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle 
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3], v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3], v_src_plu3[3], v_src_min3[3];
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	value_t a_fld, a_src;//SWG场与源的公共面面积
	value_t l_fld, l_src;//RWG场与源的公共边长度
	value_t epsp, epsm;//源处介电常数

					   //Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm_f, vm_s, vm, v_ct1, v_ct2, v_ct3;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f, 0.0f);
	Qcout << '\n' << "Fill the matrix ZDM:";
	//填充ZDM矩阵
	tool::BarAndPercent bar_perc3;
	coef = Z0 * Z0*k_ / (J0*PI4);
	for (int f = 0; f < unknowns_e; f++)
	{
		bar_perc3(f + 1, unknowns_e);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
		tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
		tri_min.getVertex(tri_min1, tri_min2, tri_min3);
		for (int s = 0; s < unknowns_t; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);
			if (vm_s != -1)
			{
				v_src_plu = mesh_ptr_->getVertex(vp);
				v_src_min = mesh_ptr_->getVertex(vm_s);

				tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);

				tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);
			}
			else
			{
				v_src_plu = mesh_ptr_->getVertex(vp);
				tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
				tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);
			}
			//fill
			if ((tri_plu1 == v_ct1 && tri_plu2 == v_ct2 && tri_plu3 == v_ct3)
				|| (tri_plu1 == v_ct1 && tri_plu2 == v_ct3 && tri_plu3 == v_ct2)
				|| (tri_plu1 == v_ct2 && tri_plu2 == v_ct1 && tri_plu3 == v_ct3)
				|| (tri_plu1 == v_ct2 && tri_plu2 == v_ct3 && tri_plu3 == v_ct1)
				|| (tri_plu1 == v_ct3 && tri_plu2 == v_ct1 && tri_plu3 == v_ct2)
				|| (tri_plu1 == v_ct3 && tri_plu2 == v_ct2 && tri_plu3 == v_ct1))
			{
				Zpp = MDZppSingular(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
				Zpp = MDZppkernel(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);


			if ((tri_min1 == v_ct1 && tri_min2 == v_ct2 && tri_min3 == v_ct3)
				|| (tri_min1 == v_ct1 && tri_min2 == v_ct3 && tri_min3 == v_ct2)
				|| (tri_min1 == v_ct2 && tri_min2 == v_ct1 && tri_min3 == v_ct3)
				|| (tri_min1 == v_ct2 && tri_min2 == v_ct3 && tri_min3 == v_ct1)
				|| (tri_min1 == v_ct3 && tri_min2 == v_ct1 && tri_min3 == v_ct2)
				|| (tri_min1 == v_ct3 && tri_min2 == v_ct2 && tri_min3 == v_ct1))
			{
				Zmp = MDZmpSingular(v_fld_min3, v_src_plu4, v_fld_min, v_fld_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
				Zmp = MDZmpkernel(v_fld_min3, v_src_plu4, v_fld_min, v_fld_plu, v_src_ct3, vm_s, epsp, epsm);

			if (vm_s != -1)
			{
				Zpm = MDZpmkernel(v_fld_plu3, v_src_min4, v_fld_plu, v_src_min, v_src_ct3, vm_s, epsp, epsm);
				Zmm = MDZmmkernel(v_fld_min3, v_src_min4, v_fld_min, v_fld_min, v_src_ct3, vm_s, epsp, epsm);
			}

			Z(f + unknowns_t, s) = coef * l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void SED_CM_FMM::fillZMM(int &start_f, int &end_f, int &start_s, int &end_s)
{
	/*int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3], v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3], v_src_plu3[3], v_src_min3[3];
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	value_t a_fld, a_src;//SWG场与源的公共面面积
	value_t l_fld, l_src;//RWG场与源的公共边长度
	value_t epsp, epsm;//源处介电常数

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm_f, vm_s, vm, v_ct1, v_ct2, v_ct3;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f, 0.0f);*/

	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	Triangle tri_plu, tri_min;
	value_t l_fld, l_src;//场与源的公共边长度

	Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm;
	Complex coef = (J0 * k_ * Z0) / PI4;
	Qcout << '\n' << "Fill the matrix ZMM:";
	//填充ZMM矩阵
	tool::BarAndPercent bar_perc4;
	//coef = (J0*k_*Z0) / PI4;
	for (int f = 0; f < unknowns_e; f++)
	{
		bar_perc4(f + 1, unknowns_e);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
		for (int s = 0; s < unknowns_e; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);

			//field face+ <-->  source face+
			if (f_fld_plu == f_src_plu)
			{
				Zpp = MMZppSingular(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			else
			{
				Zpp = MMZppKernel(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			//field face+ <--> source face-
			if (f_fld_plu == f_src_min)
			{
				Zpm = MMZpmSingular(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
			}
			else
			{
				Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
			}
			//field face- <--> source face+
			if (f_fld_min == f_src_plu)
			{
				Zmp = MMZmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
			}
			else
			{
				Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
			}
			//field face- <--> source face-
			if (f_fld_min == f_src_min)
			{
				Zmm = MMZmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
			}
			else
			{
				Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
			}
			Z(f + unknowns_t, s + unknowns_t) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void SED_CM_FMM::fillV()
{
	int nvp, nvm, fp, fm, bedge;
	VectorR3 vp, vm, vp4[4], vm4[4];
	VectorR3 vp3[3], vm3[3];
	value_t area, ln;
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	VectorR3 trans_vector = mesh_ptr_->getSize();

	tool::BarAndPercent bar_perc1;   //
	for (int u = 0; u < unknowns_t; u++)
	{
		bar_perc1(u + 1, unknowns_t);    //

		ct_ptr_->getCommonTriangle(u, nvp, nvm, fp, fm, area);
		bedge = ct_ptr_->getbedge(u);
		vp = mesh_ptr_->getVertex(nvp);
		tet_plu = mesh_ptr_->getTetrahedronRef(fp);
		tet_plu.GetVertices(vp4[0], vp4[1], vp4[2], vp4[3]);
		if (nvm != -1)
		{
			vm = mesh_ptr_->getVertex(nvm);
			tet_min = mesh_ptr_->getTetrahedronRef(fm);
			tet_min.GetVertices(vm4[0], vm4[1], vm4[2], vm4[3]);

			if (bedge == 1)
			{
				vm.x += trans_vector.x;
				vm4[0].x += trans_vector.x;
				vm4[1].x += trans_vector.x;
				vm4[2].x += trans_vector.x;
				vm4[3].x += trans_vector.x;
			}
			else if (bedge == 2)
			{
				vm.y += trans_vector.y;
				vm4[0].y += trans_vector.y;
				vm4[1].y += trans_vector.y;
				vm4[2].y += trans_vector.y;
				vm4[3].y += trans_vector.y;
			}
		}

		Complex vk = DVKernel(vp4, vm4, vp, vm, nvm);
		V(u) = area * vk / 3.0f;
	}

	tool::BarAndPercent bar_perc2;   //
	for (int u = 0; u < unknowns_e; u++)
	{
		bar_perc2(u + 1, unknowns_e);    //

		ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, ln);
		vp = mesh_ptr_->getVertex(nvp);
		vm = mesh_ptr_->getVertex(nvm);
		tri_plu = mesh_ptr_->getTriangleRef(fp);
		tri_min = mesh_ptr_->getTriangleRef(fm);
		tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
		tri_min.getVertex(vm3[0], vm3[1], vm3[2]);
		if (bedge == 1)
		{
			vm.x += trans_vector.x;
			vm3[0].x += trans_vector.x;
			vm3[1].x += trans_vector.x;
			vm3[2].x += trans_vector.x;
		}
		else if (bedge == 2)
		{
			vm.y += trans_vector.y;
			vm3[0].y += trans_vector.y;
			vm3[1].y += trans_vector.y;
			vm3[2].y += trans_vector.y;
		}

		Complex vk = MVKernel(vp3, vm3, vp, vm);
		V(u + unknowns_t) = 0.5f * ln * vk;
	}
}

void SED_CM_FMM::prepareVIE(value_t _eps1, value_t _mu1)
{
	eps1_ = _eps1;
	mu1_ = _mu1;
}

bool SED_CM_FMM::readExcEdges(const Qstring & rad_file, Qstring& stateInfo)
{
	Qifstream radStream(rad_file, std::ios::in);
	if (radStream.fail())
	{
		stateInfo.assign("read " + rad_file + " failed");
		return false;
	}
	std::pair<int, int> edge;
	while (radStream >> edge.first >> edge.second)
	{
		--edge.first;
		--edge.second;
		exc_edge_.push_back(edge);
	}

	if (radStream.eof())
		return true;
	else
		return false;
}

void SED_CM_FMM::readExcEdges(const Qstring & rad_file)
{
	Qifstream rad_stream(rad_file, std::ios::in);
	if (rad_stream.fail())
		throw component::FileError("fail loading rad file: " + rad_file);

	std::pair<int, int> edge;
	while (rad_stream >> edge.first >> edge.second)
	{
		--edge.first;
		--edge.second;
		exc_edge_.push_back(edge);
	}
	rad_stream.close();
}

void SED_CM_FMM::radiateV()
{
	int nv1, nv2;

	tool::BarAndPercent bar_perc;
	/*for (int u = 0; u < unknowns_; ++u)
	{
	bar_perc(u + 1, unknowns_);

	ce_ptr_->getCommonEdge(u, nv1, nv2);
	auto length = ce_ptr_->getCommonEdgeLength(u);

	for (const auto& elem : exc_edge_)
	{
	if ((elem.first == nv1 && elem.second == nv2) || (elem.first == nv2 && elem.second == nv1))
	{
	if (elem.first == nv1)
	V(u) = length;
	else
	V(u) = -length;
	break;
	}
	}
	}*/
}

value_t SED_CM_FMM::getBiRCS(const VectorR3& sca_k)
{
	int nvp, nvm, fp, fm;
	int nv1, nv2;
	int start_t, end_t, start_e, end_e;
	int sedid, bedge;
	VectorR3 vp, vm, vp4[4], vm4[4], vp3[3], vm3[3];
	value_t S, L, epsp = 0.0f, epsm = 0.0f;

	VectorR3 vgp_d[5], vgm_d[5];
	VectorR3 vgp_m[7], vgm_m[7];

	VectorR3 rou_p, rou_m;
	Complex G0(0, 0);
	VectorC3 Es_d(0, 0, 0), Es_m(0, 0, 0);

	for (int x = 0; x < Array_x; x++)
	{
		for (int y = 0; y < Array_y; y++)
		{
			int id = ArrayToSED(x, y);
			int blk = x * Array_y + y;
			ct_ptr_->getSEDCommonTriangleSize(id, start_t, end_t);
			ce_ptr_->getSEDCommonEdgeSize(id, start_e, end_e);
			//constructI(x, y, id);
			//auto& I_temp = I_tot[blk];
			auto& I_temp = I_tot_eigen[blk];
			for (int u = start_t; u < end_t; u++)
			{
				ct_ptr_->getCommonTriangle(u, nvp, nvm, fp, fm, S, epsp, epsm, sedid, bedge);

				vp = mesh_ptr_->getVertex(nvp);
				const auto& tet_plu = mesh_ptr_->getTetrahedronRef(fp);
				tet_plu.GetVertices(vp4[0], vp4[1], vp4[2], vp4[3]);
				shiftcoordv1(x, y, vp);
				shiftcoordv4(x, y, vp4);
				Gauss5Point(vp4[0], vp4[1], vp4[2], vp4[3], vgp_d);

				if (nvm != -1)
				{
					vm = mesh_ptr_->getVertex(nvm);
					const auto& tet_min = mesh_ptr_->getTetrahedronRef(fm);
					tet_min.GetVertices(vm4[0], vm4[1], vm4[2], vm4[3]);
					shiftcoordv1(x, y, vm);
					shiftcoordv4(x, y, vm4);
					if (bedge != -1)
					{
						shiftbedgecoordv1(bedge, vm);
						shiftbedgecoordv4(bedge, vm4);
					}
					Gauss5Point(vm4[0], vm4[1], vm4[2], vm4[3], vgm_d);
				}

				VectorC3 Esp(0, 0, 0), Esm(0, 0, 0);
				for (int g = 0; g < 5; g++)
				{
					rou_p = vgp_d[g] - vp;
					G0 = exp(J0 * k_ * (vgp_d[g] ^ sca_k));
					Esp += rou_p * wt5[g] * G0*(1.0f - (1.0f / epsp));
					if (nvm != -1)
					{
						rou_m = vm - vgm_d[g];
						G0 = exp(J0 * k_ * (vgm_d[g] ^ sca_k));
						Esm += rou_m * wt5[g] * G0*(1.0f - (1.0f / epsm));
					}

				}
				Es_d += (Esp + Esm) * (I_temp(u - start_t) * S);
			}

			for (int u = start_e; u < end_e; u++)
			{
				ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, L, sedid, bedge);
				vp = mesh_ptr_->getVertex(nvp);
				vm = mesh_ptr_->getVertex(nvm);
				const auto& tri_plu = mesh_ptr_->getTriangleRef(fp);
				const auto& tri_min = mesh_ptr_->getTriangleRef(fm);
				tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
				tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

				shiftcoordv1(x, y, vp, scale_x);
				shiftcoordv1(x, y, vm, scale_x);
				shiftcoordv3(x, y, vp3, scale_x);
				shiftcoordv3(x, y, vm3, scale_x);
				L *= pow(scale_x, x);

				if (bedge != -1)
				{
					shiftbedgecoordv1(bedge, vm);
					shiftbedgecoordv3(bedge, vm3);
				}
				Gauss7Point(vp3[0], vp3[1], vp3[2], vgp_m);
				Gauss7Point(vm3[0], vm3[1], vm3[2], vgm_m);

				VectorC3 Esp(0, 0, 0), Esm(0, 0, 0);
				for (int g = 0; g < 7; g++)
				{
					rou_p = vgp_m[g] - vp;
					G0 = exp(J0 * k_ * (vgp_m[g] ^ sca_k));
					Esp = Esp + rou_p * w7[g] * G0;

					rou_m = vm - vgm_m[g];
					G0 = exp(J0 * k_ * (vgm_m[g] ^ sca_k));
					Esm = Esm + rou_m * w7[g] * G0;
				}

				Es_m += (Esp + Esm) * (I_temp(u - start_e + end_t - start_t) * L);


			}
		}
	}



	//auto es = sca_k * (sca_k*((Z0*Es_d / 3.0f) + (Es_m / 2.0f)));
	//auto es = sca_k * (sca_k*((J0*omiga*Es_d / 3.0f) + (Es_m / 2.0f)));
	auto es = sca_k * (sca_k*((Es_d / 3.0f) + (Es_m / 2.0f)));
	value_t rcs = k_ * k_ * Z0*Z0* es.norm() / PI4;

	//return rcs;
	return 10 * log10(rcs);
}

value_t SED_CM_FMM::testRCS(const VectorR3& sca_k)
{
	int nvp, nvm, fp, fm;
	int start_t, end_t, start_e, end_e;
	int sedid, bedge;
	VectorR3 vp, vm, vp4[4], vm4[4], vp3[3], vm3[3];
	value_t S, L, epsp = 0.0f, epsm = 0.0f;

	VectorR3 vgp_d[5], vgm_d[5];
	VectorR3 vgp_m[7], vgm_m[7];

	VectorR3 rou_p, rou_m;
	Complex G0(0, 0);
	VectorC3 Es_d(0, 0, 0), Es_m(0, 0, 0);

	//Qcx_vec I_test = I_sed.col(0);
	CVector I_test = I_sed_eigen.col(0);
	for (int id = 0; id < 9; id++)
	{
		ct_ptr_->getSEDCommonTriangleSize(id, start_t, end_t);
		ce_ptr_->getSEDCommonEdgeSize(id, start_e, end_e);

		for (int u = start_t; u < end_t; u++)
		{
			ct_ptr_->getCommonTriangle(u, nvp, nvm, fp, fm, S, epsp, epsm, sedid, bedge);

			vp = mesh_ptr_->getVertex(nvp);
			const auto& tet_plu = mesh_ptr_->getTetrahedronRef(fp);
			tet_plu.GetVertices(vp4[0], vp4[1], vp4[2], vp4[3]);

			shiftcoordv1(sedid, vp);
			shiftcoordv4(sedid, vp4);

			Gauss5Point(vp4[0], vp4[1], vp4[2], vp4[3], vgp_d);

			if (nvm != -1)
			{
				vm = mesh_ptr_->getVertex(nvm);
				const auto& tet_min = mesh_ptr_->getTetrahedronRef(fm);
				tet_min.GetVertices(vm4[0], vm4[1], vm4[2], vm4[3]);

				shiftcoordv1(sedid, vm);
				shiftcoordv4(sedid, vm4);

				if (bedge != -1)
				{
					shiftbedgecoordv1(bedge, vm);
					shiftbedgecoordv4(bedge, vm4);
				}
				Gauss5Point(vm4[0], vm4[1], vm4[2], vm4[3], vgm_d);
			}

			VectorC3 Esp(0, 0, 0), Esm(0, 0, 0);
			for (int g = 0; g < 5; g++)
			{
				rou_p = vgp_d[g] - vp;
				G0 = exp(J0 * k_ * (vgp_d[g] ^ sca_k));
				Esp += rou_p * wt5[g] * G0*(1.0f - (1.0f / epsp));
				if (nvm != -1)
				{
					rou_m = vm - vgm_d[g];
					G0 = exp(J0 * k_ * (vgm_d[g] ^ sca_k));
					Esm += rou_m * wt5[g] * G0*(1.0f - (1.0f / epsm));
				}

			}
			Es_d += (Esp + Esm) * (I_test(u + start_e) * S);
		}

		for (int u = start_e; u < end_e; u++)
		{
			ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, L, sedid, bedge);
			vp = mesh_ptr_->getVertex(nvp);
			vm = mesh_ptr_->getVertex(nvm);
			const auto& tri_plu = mesh_ptr_->getTriangleRef(fp);
			const auto& tri_min = mesh_ptr_->getTriangleRef(fm);
			tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
			tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

			shiftcoordv1(sedid, vp);
			shiftcoordv1(sedid, vm);
			shiftcoordv3(sedid, vp3);
			shiftcoordv3(sedid, vm3);

			if (bedge != -1)
			{
				shiftbedgecoordv1(bedge, vm);
				shiftbedgecoordv3(bedge, vm3);
			}
			Gauss7Point(vp3[0], vp3[1], vp3[2], vgp_m);
			Gauss7Point(vm3[0], vm3[1], vm3[2], vgm_m);

			VectorC3 Esp(0, 0, 0), Esm(0, 0, 0);
			for (int g = 0; g < 7; g++)
			{
				rou_p = vgp_m[g] - vp;
				G0 = exp(J0 * k_ * (vgp_m[g] ^ sca_k));
				Esp = Esp + rou_p * w7[g] * G0;

				rou_m = vm - vgm_m[g];
				G0 = exp(J0 * k_ * (vgm_m[g] ^ sca_k));
				Esm = Esm + rou_m * w7[g] * G0;
			}
			Es_m += (Esp + Esm) * (I_test(u + end_t) * L);
		}
	}
	auto es = sca_k * (sca_k*((Es_d / 3.0f) + (Es_m / 2.0f)));
	//auto es = sca_k * (sca_k*((Z0*Es_d / 3.0f) + (Es_m / 2.0f)));
	value_t rcs = k_ * k_ * Z0*Z0* es.norm() / PI4;

	//return rcs;
	return 10 * log10(rcs);
}

void SED_CM_FMM::getEFiled(component::FieldData * pdata, const VectorR3 & rad_k, const VectorR3 & rad_ev, const VectorR3 & rad_eh) const
{
	int nvp, nvm, fp, fm;
	int nv1, nv2;
	int start_t, end_t, start_e, end_e;
	int sedid, bedge;
	VectorR3 vp, vm, vp4[4], vm4[4], vp3[3], vm3[3];
	value_t S, L, epsp = 0.0f, epsm = 0.0f;

	VectorR3 vgp_d[5], vgm_d[5];
	VectorR3 vgp_m[7], vgm_m[7];

	VectorR3 rou_p, rou_m;
	Complex G0(0, 0);
	VectorC3 Es_d(0, 0, 0), Es_m(0, 0, 0);

	for (int x = 0; x < Array_x; x++)
	{
		for (int y = 0; y < Array_y; y++)
		{
			int id = ArrayToSED(x, y);
			int blk = x * Array_y + y;
			ct_ptr_->getSEDCommonTriangleSize(id, start_t, end_t);
			ce_ptr_->getSEDCommonEdgeSize(id, start_e, end_e);
			//constructI(x, y, id);
			//auto& I_temp = I_tot[blk];
			auto& I_temp = I_tot_eigen[blk];
			for (int u = start_t; u < end_t; u++)
			{
				ct_ptr_->getCommonTriangle(u, nvp, nvm, fp, fm, S, epsp, epsm, sedid, bedge);

				vp = mesh_ptr_->getVertex(nvp);
				const auto& tet_plu = mesh_ptr_->getTetrahedronRef(fp);
				tet_plu.GetVertices(vp4[0], vp4[1], vp4[2], vp4[3]);
				shiftcoordv1(x, y, vp);
				shiftcoordv4(x, y, vp4);
				Gauss5Point(vp4[0], vp4[1], vp4[2], vp4[3], vgp_d);

				if (nvm != -1)
				{
					vm = mesh_ptr_->getVertex(nvm);
					const auto& tet_min = mesh_ptr_->getTetrahedronRef(fm);
					tet_min.GetVertices(vm4[0], vm4[1], vm4[2], vm4[3]);
					shiftcoordv1(x, y, vm);
					shiftcoordv4(x, y, vm4);
					if (bedge != -1)
					{
						shiftbedgecoordv1(bedge, vm);
						shiftbedgecoordv4(bedge, vm4);
					}
					Gauss5Point(vm4[0], vm4[1], vm4[2], vm4[3], vgm_d);
				}

				VectorC3 Esp(0, 0, 0), Esm(0, 0, 0);
				for (int g = 0; g < 5; g++)
				{
					rou_p = vgp_d[g] - vp;
					G0 = exp(J0 * k_ * (vgp_d[g] ^ rad_k));
					Esp += rou_p * wt5[g] * G0*(1.0f - (1.0f / epsp));
					if (nvm != -1)
					{
						rou_m = vm - vgm_d[g];
						G0 = exp(J0 * k_ * (vgm_d[g] ^ rad_k));
						Esm += rou_m * wt5[g] * G0*(1.0f - (1.0f / epsm));
					}

				}
				Es_d += (Esp + Esm) * (I_temp(u - start_t) * S);
			}

			for (int u = start_e; u < end_e; u++)
			{
				ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, L, sedid, bedge);
				vp = mesh_ptr_->getVertex(nvp);
				vm = mesh_ptr_->getVertex(nvm);
				const auto& tri_plu = mesh_ptr_->getTriangleRef(fp);
				const auto& tri_min = mesh_ptr_->getTriangleRef(fm);
				tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
				tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

				shiftcoordv1(x, y, vp);
				shiftcoordv1(x, y, vm);
				shiftcoordv3(x, y, vp3);
				shiftcoordv3(x, y, vm3);

				if (bedge != -1)
				{
					shiftbedgecoordv1(bedge, vm);
					shiftbedgecoordv3(bedge, vm3);
				}
				Gauss7Point(vp3[0], vp3[1], vp3[2], vgp_m);
				Gauss7Point(vm3[0], vm3[1], vm3[2], vgm_m);

				VectorC3 Esp(0, 0, 0), Esm(0, 0, 0);
				for (int g = 0; g < 7; g++)
				{
					rou_p = vgp_m[g] - vp;
					G0 = exp(J0 * k_ * (vgp_m[g] ^ rad_k));
					Esp = Esp + rou_p * w7[g] * G0;

					rou_m = vm - vgm_m[g];
					G0 = exp(J0 * k_ * (vgm_m[g] ^ rad_k));
					Esm = Esm + rou_m * w7[g] * G0;
				}

				Es_m += (Esp + Esm) * (I_temp(u - start_e + end_t - start_t) * L);


			}
		}
	}



	//auto es = sca_k * (sca_k*((Z0*Es_d / 3.0f) + (Es_m / 2.0f)));
	//auto es = sca_k * (sca_k*((J0*omiga*Es_d / 3.0f) + (Es_m / 2.0f)));
	auto es = (Es_d / 3.0f) + (Es_m / 2.0f);
	const auto coeff = -J0 * Z0*k_ / (PI4);
	pdata->ftheta = coeff * (es^rad_ev);
	pdata->fphi = coeff * (es^rad_eh);

}

void SED_CM_FMM::getCurrentOnFeed()
{
	int nv1, nv2;
	int nvp, nvm, fp, fm;
	int start_t, end_t, start_e, end_e;
	int sedid, bedge;
	value_t S, L;

	Qofstream outputI(dir_ + "/I_feed.txt");
	for (int x = 0; x < Array_x; x++)
	{
		for (int y = 0; y < Array_y; y++)
		{
			int id = ArrayToSED(x, y);
			int blk = x * Array_y + y;

			ce_ptr_->getSEDCommonEdgeSize(id, start_e, end_e);

			auto& I_temp = I_tot_eigen[blk];
			for (int u = start_e; u < end_e; u++)
			{
				ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, L, sedid, bedge);
				ce_ptr_->getCommonEdge(u, nv1, nv2);

				for (const auto& elem : exc_edge_)
				{
					if ((elem.first == nv1 && elem.second == nv2) || (elem.first == nv2 && elem.second == nv1))
					{
						if (elem.first == nv1)
							I_feed(blk) = (I_temp(u - start_e) * L);
						else
							I_feed(blk) = (-I_temp(u - start_e) * L);
						outputI << std::setw(18) << I_feed(blk).real() << std::setw(18) << I_feed(blk).imag() << '\n';
					}
				}
			}
		}
	}

	outputI.flush();
	outputI.close();

}

void SED_CM_FMM::getNearEField(std::vector<component::NearFieldData>* data) const
{
	VectorR3 ori = nearfield_.origin_point;
	VectorR3 end = nearfield_.end_point;
	int samp_x_num = nearfield_.sampling_x;
	int samp_y_num = nearfield_.sampling_y;
	int samp_z_num = nearfield_.sampling_z;

	value_t delta_x, delta_y, delta_z;
	std::vector<VectorR3> pointArr;

	if ((end - ori).x < 1e-6)
		delta_x = 0;
	else if (samp_x_num > 1)
		delta_x = (end - ori).x / (samp_x_num - 1);
	else
		delta_x = 0;

	if ((end - ori).y < 1e-6)
		delta_y = 0;
	else if (samp_y_num > 1)
		delta_y = (end - ori).y / (samp_y_num - 1);
	else
		delta_y = 0;

	if ((end - ori).z < 1e-6)
		delta_z = 0;
	else if (samp_z_num > 1)
		delta_z = (end - ori).z / (samp_z_num - 1);
	else
		delta_z = 0;

	VectorR3 point = nearfield_.origin_point;
	VectorR3 point_temp;
	for (int i = 0; i < samp_z_num; i++)
	{
		point_temp.z = point.z + i * delta_z;
		for (int j = 0; j < samp_y_num; j++)
		{
			point_temp.y = point.y + j * delta_y;
			for (int k = 0; k < samp_x_num; k++)
			{
				point_temp.x = point.x + k * delta_x;
				pointArr.push_back(point_temp);
			}
		}
	}
	int point_num = pointArr.size();

	VectorC3 E_near;
	Complex coff = -J0 * k_*Z0 / PI4;

	/*for (int p = 0; p < point_num; p++)
	{
	VectorC3 E_vol, E_surf, E_inc, e_vol_temp, e_surf_temp;
	VectorR3 ob_point = pointArr[p];
	for (int t = 0; t < unknowns_t; t++)
	{
	NEFkernel_vol(ob_point, t, e_vol_temp);
	E_vol += I(t)*e_vol_temp;
	}
	for (int e = 0; e < unknowns_e; e++)
	{
	NEFkernel_surf(ob_point, e, e_surf_temp);
	E_surf += I(e + unknowns_t)*e_surf_temp;
	}
	E_near = coff * (E_vol + E_surf);
	E_inc = inc_e_ * exp(-J0 * k_*(ob_point^inc_k_));
	NearFieldData EFdata(ob_point, E_near);
	EFdata.sca_inc = E_near + E_inc;
	data->push_back(EFdata);
	}*/

	for (int p = 0; p < point_num; p++)
	{
		int start_t, end_t, start_e, end_e;
		VectorC3 E_vol, E_surf, E_inc, e_vol_temp, e_surf_temp;
		VectorR3 ob_point = pointArr[p];
		for (int x = 0; x < Array_x; x++)
		{
			for (int y = 0; y < Array_y; y++)
			{
				int id = ArrayToSED(x, y);
				int blk = x * Array_y + y;
				//ct_ptr_->getSEDCommonTriangleSize(id, start_t, end_t);
				ce_ptr_->getSEDCommonEdgeSize(id, start_e, end_e);
				//constructI(x, y, id);
				//auto& I_temp = I_tot[blk];
				auto& I_temp = I_tot_eigen[blk];

				for (int u = start_e; u < end_e; u++)
				{
					NEFkernel_surf(ob_point, x, y, u, e_surf_temp);
					E_surf += I_temp(u - start_e)*e_surf_temp;
				}
			}
		}
		E_near = coff * (E_surf);
		E_inc = inc_e_ * exp(-J0 * k_*(ob_point^inc_k_));
		NearFieldData EFdata(ob_point, E_near);
		EFdata.sca_inc = E_near + E_inc;
		data->push_back(EFdata);
	}

}

void SED_CM_FMM::NEFkernel_vol(VectorR3& ob, int& unk, VectorC3& e_vol) const
{
	int nvp, nvm, tp, tm, v1, v2, v3;
	value_t area, epsp, epsm;
	VectorR3 vp, vm, vp4[4], vm4[4], vt3[3];
	VectorR3 vgp4[4], vgm4[4], vg3[3];

	ct_ptr_->getCommonTriangle(unk, nvp, nvm, tp, tm, area, epsp, epsm);
	ct_ptr_->getCommonTriangle(unk, v1, v2, v3);

	auto& tet_plu = mesh_ptr_->getTetrahedronRef(tp);
	tet_plu.GetVertices(vp4[0], vp4[1], vp4[2], vp4[3]);
	Gauss4Point(vp4[0], vp4[1], vp4[2], vp4[3], vgp4);
	vp = mesh_ptr_->getVertex(nvp);

	if (tm != -1)
	{
		auto& tet_min = mesh_ptr_->getTetrahedronRef(tm);
		tet_min.GetVertices(vm4[0], vm4[1], vm4[2], vm4[3]);
		Gauss4Point(vm4[0], vm4[1], vm4[2], vm4[3], vgm4);
		vm = mesh_ptr_->getVertex(nvm);
	}

	vt3[0] = mesh_ptr_->getVertex(v1);
	vt3[1] = mesh_ptr_->getVertex(v2);
	vt3[2] = mesh_ptr_->getVertex(v3);
	Gauss3Point(vt3[0], vt3[1], vt3[2], vg3);

	VectorC3 E1, E2, E3;
	VectorR3 rou_p, rou_m, r;
	Complex G0, G_grad;
	value_t kp, km;
	kp = 1.0f - (1.0f / epsp);
	km = 1.0f - (1.0f / epsm);

	for (int i = 0; i < 4; i++)
	{
		rou_p = vgp4[i] - vp;
		r = ob - vgp4[i];

		G0 = (exp(-J0 * k_*r.Norm())) / r.Norm();
		G_grad = G0 * (-J0 * k_ - (1.0f / r.Norm()));

		E1 += wt4[i] * rou_p*kp*G0;
		E2 += wt4[i] * kp*r.Normalize()*G_grad;
	}

	if (abs(epsp - epsm) > 1e-5)
	{
		for (int i = 0; i < 3; i++)
		{
			r = ob - vg3[i];
			G0 = (exp(-J0 * k_*r.Norm())) / r.Norm();
			G_grad = G0 * (-J0 * k_ - (1.0f / r.Norm()));
			E3 += w3[i] * (km - kp) *r.Normalize()*G_grad;
		}
	}
	if (tm != -1)
	{
		for (int i = 0; i < 4; i++)
		{
			rou_m = vm - vgm4[i];
			r = vgm4[i] - ob;

			G0 = (exp(-J0 * k_*r.Norm())) / r.Norm();
			G_grad = G0 * (-J0 * k_ - (1.0f / r.Norm()));

			E1 += wt4[i] * rou_m * km * G0;
			E2 += wt4[i] * km * r.Normalize() * G_grad;
		}
	}

	e_vol = (area*((E1 / 3.0f) + (E2 + E3) / (k_*k_)));
}

void SED_CM_FMM::NEFkernel_surf(VectorR3& ob, int& unk, VectorC3& e_surf) const
{
	int nvp, nvm, fp, fm, v1, v2;
	value_t len;
	VectorR3 vp, vm, vp3[3], vm3[3];
	VectorR3 vgp3[3], vgm3[3];

	ce_ptr_->getCommonEdge(unk, nvp, nvm, fp, fm, len);

	auto &tri_plu = mesh_ptr_->getTriangleRef(fp);
	auto &tri_min = mesh_ptr_->getTriangleRef(fm);

	tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
	tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

	vp = mesh_ptr_->getVertex(nvp);
	vm = mesh_ptr_->getVertex(nvm);

	Gauss3Point(vp3[0], vp3[1], vp3[2], vgp3);
	Gauss3Point(vm3[0], vm3[1], vm3[2], vgm3);

	VectorR3 r, rou_p, rou_m;
	VectorC3 E1, E2;
	Complex G0, G_grad;

	for (int i = 0; i < 3; i++)
	{
		r = ob - vgp3[i];
		rou_p = vgp3[i] - vp;

		G0 = (exp(-J0 * k_*r.Norm())) / r.Norm();
		G_grad = G0 * (-J0 * k_ - (1.0f / r.Norm()));

		E1 += w3[i] * rou_p*G0;
		E2 += w3[i] * r.Normalize()*G_grad;
	}

	for (int i = 0; i < 3; i++)
	{
		r = vgm3[i] - ob;
		rou_m = vm - vgm3[i];

		G0 = (exp(-J0 * k_*r.Norm())) / r.Norm();
		G_grad = G0 * (-J0 * k_ - (1.0f / r.Norm()));

		E1 += w3[i] * rou_m*G0;
		E2 += w3[i] * r.Normalize() *G_grad;
	}

	e_surf = (len*((E1 / 2.0f) + E2 / (k_*k_)));
}

void SED_CM_FMM::NEFkernel_surf(VectorR3& ob, int& _x, int& _y, int& unk, component::VectorC3& e_surf) const
{
	int nvp, nvm, fp, fm, v1, v2, sedid, bedge;
	value_t L;
	VectorR3 vp, vm, vp3[3], vm3[3];
	VectorR3 vgp3[3], vgm3[3];


	VectorR3 r, rou_p, rou_m;
	VectorC3 E1, E2;
	Complex G0, G_grad;
	value_t scale = scale_x;
	ce_ptr_->getCommonEdge(unk, nvp, nvm, fp, fm, L, sedid, bedge);
	vp = mesh_ptr_->getVertex(nvp);
	vm = mesh_ptr_->getVertex(nvm);
	const auto& tri_plu = mesh_ptr_->getTriangleRef(fp);
	const auto& tri_min = mesh_ptr_->getTriangleRef(fm);
	tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
	tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

	shiftcoordv1(_x, _y, vp, scale);
	shiftcoordv1(_x, _y, vm, scale);
	shiftcoordv3(_x, _y, vp3, scale);
	shiftcoordv3(_x, _y, vm3, scale);
	L *= pow(scale, _x);
	if (bedge != -1)
	{
		shiftbedgecoordv1(bedge, vm);
		shiftbedgecoordv3(bedge, vm3);
	}
	Gauss3Point(vp3[0], vp3[1], vp3[2], vgp3);
	Gauss3Point(vm3[0], vm3[1], vm3[2], vgm3);

	for (int i = 0; i < 3; i++)
	{
		r = ob - vgp3[i];
		rou_p = vgp3[i] - vp;

		G0 = (exp(-J0 * k_*r.Norm())) / r.Norm();
		G_grad = G0 * (-J0 * k_ - (1.0f / r.Norm()));

		E1 += w3[i] * rou_p*G0;
		E2 += w3[i] * r.Normalize()*G_grad;
	}

	for (int i = 0; i < 3; i++)
	{
		r = vgm3[i] - ob;
		rou_m = vm - vgm3[i];

		G0 = (exp(-J0 * k_*r.Norm())) / r.Norm();
		G_grad = G0 * (-J0 * k_ - (1.0f / r.Norm()));

		E1 += w3[i] * rou_m*G0;
		E2 += w3[i] * r.Normalize() *G_grad;
	}

	e_surf = (L*((E1 / 2.0f) + E2 / (k_*k_)));
}

value_t SED_CM_FMM::getReflectField(const VectorR3& sca_k) const
{
	int nvp, nvm, fp, fm;
	VectorR3 vp, vm, vp4[4], vm4[4], vp3[3], vm3[3];
	value_t S, L, epsp = 0.0f, epsm = 0.0f;

	VectorR3 vgp_d[5], vgm_d[5];
	VectorR3 vgp_m[7], vgm_m[7];

	VectorR3 rou_p, rou_m;
	Complex G0(0, 0);
	VectorC3 Es_d(0, 0, 0), Es_m(0, 0, 0);

	for (int u = 0; u < unknowns_t; u++)
	{
		ct_ptr_->getCommonTriangle(u, nvp, nvm, fp, fm, S, epsp, epsm);

		vp = mesh_ptr_->getVertex(nvp);
		const auto& tet_plu = mesh_ptr_->getTetrahedronRef(fp);
		tet_plu.GetVertices(vp4[0], vp4[1], vp4[2], vp4[3]);
		Gauss5Point(vp4[0], vp4[1], vp4[2], vp4[3], vgp_d);

		if (nvm != -1)
		{
			vm = mesh_ptr_->getVertex(nvm);
			const auto& tet_min = mesh_ptr_->getTetrahedronRef(fm);
			tet_min.GetVertices(vm4[0], vm4[1], vm4[2], vm4[3]);
			Gauss5Point(vm4[0], vm4[1], vm4[2], vm4[3], vgm_d);
		}

		VectorC3 Esp(0, 0, 0), Esm(0, 0, 0);
		for (int g = 0; g < 5; g++)
		{
			rou_p = vgp_d[g] - vp;
			G0 = exp(J0 * k_ * (vgp_d[g] ^ sca_k));
			Esp += rou_p * wt5[g] * G0*(1.0f - (1.0f / epsp));
			if (nvm != -1)
			{
				rou_m = vm - vgm_d[g];
				G0 = exp(J0 * k_ * (vgm_d[g] ^ sca_k));
				Esm += rou_m * wt5[g] * G0*(1.0f - (1.0f / epsm));
			}

		}
		Es_d += (Esp + Esm) * (I[u] * S);
	}

	for (int u = 0; u < unknowns_e; u++)
	{
		ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, L);
		vp = mesh_ptr_->getVertex(nvp);
		vm = mesh_ptr_->getVertex(nvm);
		const auto& tri_plu = mesh_ptr_->getTriangleRef(fp);
		const auto& tri_min = mesh_ptr_->getTriangleRef(fm);
		tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
		tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

		Gauss7Point(vp3[0], vp3[1], vp3[2], vgp_m);
		Gauss7Point(vm3[0], vm3[1], vm3[2], vgm_m);

		VectorC3 Esp(0, 0, 0), Esm(0, 0, 0);
		for (int g = 0; g < 7; g++)
		{
			rou_p = vgp_m[g] - vp;
			G0 = exp(J0 * k_ * (vgp_m[g] ^ sca_k));
			Esp = Esp + rou_p * w7[g] * G0;

			rou_m = vm - vgm_m[g];
			G0 = exp(J0 * k_ * (vgm_m[g] ^ sca_k));
			Esm = Esm + rou_m * w7[g] * G0;
		}
		Es_m += (Esp + Esm) * (I[u + unknowns_t] * L);
	}


	auto es = sca_k * (sca_k*((Z0*Es_d / 3.0f) + (Es_m / 2.0f)));
	value_t rcs = k_ * k_ * Z0*Z0* es.norm() / PI4;

	return 10 * log10(rcs);
}

bool SED_CM_FMM::writeZIVData()
{

	Qofstream outputZ(dir_ + "/martix_Z_RED.txt");
	outputZ << Z_red_eigen;
	outputZ.flush();
	outputZ.close();

	Qofstream outputV_red(dir_ + "/martix_V_RED.txt");
	outputV_red << V_red_eigen;
	outputV_red.flush();
	outputV_red.close();

	Qofstream outputI_first(dir_ + "/martix_I_first.txt");
	outputI_first << I_tot_eigen[0];
	outputI_first.flush();
	outputI_first.close();

	Qofstream outputV_sed(dir_ + "/martix_V_SED.txt");
	outputV_sed << V_sed_eigen;
	outputV_sed.flush();
	outputV_sed.close();

	Qofstream outputnearZ_(dir_ + "/martix_nearZ_.txt");
	outputnearZ_ << nearZ_;
	outputnearZ_.flush();
	outputnearZ_.close();
	return true;
}

void SED_CM_FMM::calculateSurfaceCurrent(std::vector<component::CurrentData>* currents) const
{

	int tri_num, start, end;
	int nv1, nv2, nvp, nvm, fp, fm;
	value_t len;
	VectorR3 vf[3], vs[3];
	tri_num = mesh_ptr_->getTriangleNum();
	currents->resize(Array_x*Array_y*tri_num);
	value_t scale = scale_x;

	for (int x = 0; x < Array_x; x++)
	{
		for (int y = 0; y < Array_y; y++)
		{
			int bedge;
			int sed_id = ArrayToSED(x, y);
			int ver_num = mesh_ptr_->getNodeNum();

			std::vector<std::vector<value_t>> cur_ver(ver_num);
			std::vector<std::vector<int>> tri_rwg(tri_num);
			std::vector<value_t> cur_ver_mag(ver_num);
			std::vector<value_t> cur_cen(tri_num);

			ce_ptr_->getSEDCommonEdgeSize(sed_id, start, end);
			CVector I_vec = I_tot_eigen[y + x * Array_y];

			//遍历公共边，将公共边压入至三角形中，之后可知每个三角形的公共边
			for (int u = start; u < end; u++)
			{
				ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, len, bedge);

				tri_rwg[fp].push_back(u);
				if (bedge == -1)
					tri_rwg[fm].push_back(u);
			}

			//遍历三角形，将每个三角形上的所有公共边对三个顶点以及中心点的贡献求和压入顶点电流栈
			for (int t = 0; t < tri_num; t++)
			{
				int n[3];
				VectorR3 v[3];

				Triangle tri = mesh_ptr_->getTriangleRef(t);
				tri.getVertex(n[0], n[1], n[2]);
				v[0] = mesh_ptr_->getVertex(n[0]);
				v[1] = mesh_ptr_->getVertex(n[1]);
				v[2] = mesh_ptr_->getVertex(n[2]);
				scalecoordv3(x, y, v, scale);

				value_t area = Area(v[0], v[1], v[2]);
				VectorR3 v_cen = (v[0] + v[1] + v[2]) / 3.0f;

				auto& ce = tri_rwg[t];

				for (int k = 0; k < 3; k++)//遍历该三角形的三个顶点
				{
					VectorC3 cur_node_vec(0.0f, 0.0f, 0.0f);
					//VectorC3 cur_cen_vec(0.0f, 0.0f, 0.0f);
					for (int e = 0; e < ce.size(); e++)//遍历该三角形的所有公共边，计算公共边对该顶点的贡献
					{
						ce_ptr_->getCommonEdge(ce[e], nv1, nv2);
						ce_ptr_->getCommonEdge(ce[e], nvp, nvm, fp, fm, len);
						len *= pow(scale, x);

						value_t sign = 1.0f;
						VectorR3 vx;

						if (t == fp)
						{
							vx = mesh_ptr_->getVertex(nvp);
							scalecoordv1(x, y, vx, scale);
						}
						else
						{
							vx = mesh_ptr_->getVertex(nvm);
							scalecoordv1(x, y, vx, scale);
							sign = -1.0f;
						}

						if (n[k] == nv1)
						{
							VectorR3 v1 = mesh_ptr_->getVertex(n[k]);
							scalecoordv1(x, y, v1, scale);
							cur_node_vec += sign * I_vec(ce[e] - start)*(v1 - vx)*len / (2.0f * area);
						}
						else if (n[k] == nv2)
						{
							VectorR3 v2 = mesh_ptr_->getVertex(n[k]);
							scalecoordv1(x, y, v2, scale);
							cur_node_vec += sign * I_vec(ce[e] - start)*(v2 - vx)*len / (2.0f * area);
						}
					}
					cur_ver[n[k]].push_back(std::sqrt(cur_node_vec.norm()));
				}

				VectorC3 cur_cen_vec(0.0f, 0.0f, 0.0f);
				for (int e = 0; e < ce.size(); e++)//计算所有公共边对中心点的贡献
				{
					ce_ptr_->getCommonEdge(ce[e], nv1, nv2);
					ce_ptr_->getCommonEdge(ce[e], nvp, nvp, fp, fm, len);
					len *= pow(scale, x);

					if (t == fp)
					{
						VectorR3 vp = mesh_ptr_->getVertex(nvp);
						scalecoordv1(x, y, vp, scale);
						cur_cen_vec += I_vec(ce[e] - start)*(v_cen - vp)*len / (2.0f * area);
					}
					else
					{
						VectorR3 vm = mesh_ptr_->getVertex(nvm);
						scalecoordv1(x, y, vm, scale);
						cur_cen_vec += I_vec(ce[e] - start)*(vm - v_cen)*len / (2.0f * area);
					}
				}
				cur_cen[t] = (std::sqrt(cur_cen_vec.norm()));
			}

			//遍历所有顶点，将顶点的电流值求和取平均
			for (int n = 0; n < ver_num; n++)
			{
				cur_ver_mag[n] = cur_ver[n].empty() ? 0.0f : std::accumulate(cur_ver[n].begin(), cur_ver[n].end(), 0.0f) / cur_ver[n].size();
			}

			//遍历三角形，将三角形的相关信息以及电流信息写入文件中
			for (int t = 0; t < tri_num; t++)
			{
				auto& data = (*currents)[t + (x*Array_y + y)*tri_num];
				auto& tri = mesh_ptr_->getTriangleRef(t);
				int n[3];
				VectorR3 v[3];
				tri.getVertex(n[0], n[1], n[2]);
				v[0] = mesh_ptr_->getVertex(n[0]);
				v[1] = mesh_ptr_->getVertex(n[1]);
				v[2] = mesh_ptr_->getVertex(n[2]);

				//value_t scale = scale_x;
				shiftcoordv3(x, y, v, scale);
				data.magn1 = cur_ver_mag[n[0]];
				data.magn2 = cur_ver_mag[n[1]];
				data.magn3 = cur_ver_mag[n[2]];
				data.magnc = cur_cen[t];
				data.v1 = v[0];
				data.v2 = v[1];
				data.v3 = v[2];
			}

			cur_ver.clear();
			tri_rwg.clear();
			cur_ver_mag.clear();
			cur_cen.clear();
		}
	}
}

void SED_CM_FMM::readFEKOcurrent()
{
	int tri_num;
	auto maxSize = std::numeric_limits<std::streamsize>::max();
	Qifstream input_stream(FEKOcur_path, std::ios::in);

	//if (input_stream.fail())
	//return false;

	for (int t = 0; t < 9; ++t)
		input_stream.ignore(maxSize, '\n');

	input_stream.ignore(maxSize, ':');
	input_stream >> tri_num;
	input_stream.get();
	input_stream.ignore(maxSize, '\n');
	input_stream.ignore(maxSize, '\n');
	for (int i = 0; i < tri_num; i++)
	{
		VectorR3 cencoord;
		Complex jcen;
		value_t jx_re, jx_im, jy_re, jy_im, jz_re, jz_im;
		int id;
		input_stream >> id >> cencoord.x >> cencoord.y >> cencoord.z >> jx_re >> jx_im >> jy_re >> jy_im >> jz_re >> jz_im;
		input_stream.ignore(maxSize, '\n');
		Complex jx(jx_re, jx_im);
		Complex jy(jy_re, jy_im);
		Complex jz(jz_re, jz_im);
		jcen = jx + jy + jz;
		Centercoord.push_back(cencoord);
		J_center.push_back(jcen);
	}
}

void SED_CM_FMM::fillSEDZ()
{
	//int temp_t = 0, temp_e = 0;
	//int bedge_t, sedid_t, bedge_e, sedid_e;
	tool::BarAndPercent bar_perc1;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;
	int p = 0;

	//tool::BarAndPercent bar_perc1;
	for (int f_sed = 0; f_sed < 9; f_sed++)
	{
		ct_ptr_->getSEDCommonTriangleSize(f_sed, start_tf, end_tf);
		ce_ptr_->getSEDCommonEdgeSize(f_sed, start_ef, end_ef);
		for (int s_sed = 0; s_sed < 9; s_sed++)
		{
			ct_ptr_->getSEDCommonTriangleSize(s_sed, start_ts, end_ts);
			ce_ptr_->getSEDCommonEdgeSize(s_sed, start_es, end_es);
			fillSEDZDD(f_sed, s_sed, start_tf, end_tf, start_ts, end_ts, start_ef, end_ef, start_es, end_es);
			fillSEDZDM(f_sed, s_sed, start_tf, end_tf, start_ts, end_ts, start_ef, end_ef, start_es, end_es);
			fillSEDZMD(f_sed, s_sed, start_tf, end_tf, start_ts, end_ts, start_ef, end_ef, start_es, end_es);
			fillSEDZMM(f_sed, s_sed, start_tf, end_tf, start_ts, end_ts, start_ef, end_ef, start_es, end_es);
			bar_perc1(p + 1, 81);
			p = p + 1;
		}
	}
}

void SED_CM_FMM::fillSEDZDD_MP(int &_f_sed, int &_s_sed, int &_start_tf, int &_end_tf, int &_start_ts, int &_end_ts, int &_start_ef, int &_end_ef, int &_start_es, int &_end_es)
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3];
	Tetrahedron tet_plu, tet_min;
	value_t a_fld, a_src;//场与源的公共面面积
	value_t epsp, epsm;//源处介电常数
	int sedid_f, bedge_f, sedid_s, bedge_s;

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp_f, vp_s, vm_f, vm_s, v_ct1, v_ct2, v_ct3;

	Complex coef(0.0f, 0.0f);
	KIM ik = SED_KernelIntegral(_f_sed, _s_sed);
	//Qcout << '\n' << "Fill the matrix ZDD:";
	//填充ZDD矩阵
	//tool::BarAndPercent bar_perc1;
	//coef = Z0 * Z0 / (J0* PI4);
	//coef = Z0 * omiga / (PI4);
	coef = Z0 / (J0* PI4);
	//#pragma omp parallel
	for (int f = _start_tf; f < _end_tf; f++)
	{
		//bar_perc1(f-_start_tf + 1, _end_tf-_start_tf);

		ct_ptr_->getCommonTriangle(f, vp_f, vm_f, t_fld_plu, t_fld_min, a_fld, sedid_f, bedge_f);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);

		shiftcoordv3(sedid_f, v_fld_ct3);//平移公共面三角形坐标

		v_fld_plu = mesh_ptr_->getVertex(vp_f);//场点+
		tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
		tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);

		shiftcoordv1(sedid_f, v_fld_plu);//平移正顶点坐标
		shiftcoordv4(sedid_f, v_fld_plu4);//平移正四面体四点坐标

		int x1 = sedid_f / 3;
		int y1 = sedid_f % 3;

		if (vm_f != -1)
		{
			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);


			shiftcoordv1(sedid_f, v_fld_min);
			shiftcoordv4(sedid_f, v_fld_min4);

			if (bedge_f != -1)//多一次平移边界处四面体
			{
				shiftbedgecoordv1(bedge_f, v_fld_min);
				shiftbedgecoordv4(bedge_f, v_fld_min4);
			}
		}

		for (int s = _start_ts; s < _end_ts; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp_s, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm, sedid_s, bedge_s);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);

			shiftcoordv3(sedid_s, v_src_ct3);

			v_src_plu = mesh_ptr_->getVertex(vp_s);
			tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
			tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);

			shiftcoordv1(sedid_s, v_src_plu);
			shiftcoordv4(sedid_s, v_src_plu4);

			int x2 = sedid_s / 3;
			int y2 = sedid_s % 3;

			if (vm_s != -1)
			{
				v_src_min = mesh_ptr_->getVertex(vm_s);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);

				shiftcoordv1(sedid_s, v_src_min);
				shiftcoordv4(sedid_s, v_src_min4);

				if (bedge_s != -1)
				{
					shiftbedgecoordv1(bedge_s, v_src_min);
					shiftbedgecoordv4(bedge_s, v_src_min4);
				}
			}
			if (bedge_f == -1 && bedge_s == -1)
			{
				if (sedid_f == sedid_s)
				{
					if (t_fld_plu == t_src_plu)
						Zpp = DDZppSingular_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
					else
						Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_plu == t_src_min)
						Zpm = DDZpmSingular_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);
					else if (t_src_min != -1)
						Zpm = DDZpmKernel_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_min == t_src_plu)
						Zmp = DDZmpSingular_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);
					else if (t_fld_min != -1)
						Zmp = DDZmpKernel_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_min != -1 && t_src_min != -1)
					{
						if (t_fld_min == t_src_min)
							Zmm = DDZmmSingular_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
						else
							Zmm = DDZmmKernel_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
					}

				}
				else
				{
					Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);

					if (t_src_min != -1)
						Zpm = DDZpmKernel_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);

					if (t_fld_min != -1)
						Zmp = DDZmpKernel_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);

					if (t_fld_min != -1 && t_src_min != -1)
						Zmm = DDZmmKernel_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
				}
			}
			else
			{
				//field tetrahedron+ <-->source tetrahedron+ 
				if ((t_fld_plu == t_src_plu) && (sedid_f == sedid_s))
				{
					Zpp = DDZppSingular_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
				}
				else
				{
					Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
				}
				//field tetrahedron+ <-->source tetrahedron-
				if (t_fld_plu == t_src_min)
				{
					if ((centerv4(v_fld_plu4) - centerv4(v_src_min4)).Norm() < 1e-5)
					{
						Zpm = DDZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
					}
					else
					{
						Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
					}
				}
				else if (t_src_min != -1)
				{
					Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
				}
				//field tetrahedron- <-->source tetrahedron+
				if (t_fld_min == t_src_plu)
				{
					if ((centerv4(v_fld_min4) - centerv4(v_src_plu4)).Norm() < 1e-5)
					{
						Zmp = DDZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
					}
					else
					{
						Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
					}
				}
				else if (t_fld_min != -1)
				{
					Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
				}
				//field tetrahedron- <-->source tetrahedron-
				if ((t_fld_min != -1) && (t_src_min != -1))
				{
					if (t_fld_min == t_src_min)
					{
						if ((centerv4(v_fld_min4) - centerv4(v_src_min4)).Norm() < 1e-5)
						{
							Zmm = DDZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
						}
						else
						{
							Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
						}
					}
					else
						Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
				}
			}


			Z_sed(_start_ef + f, _start_es + s) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
			Z_sed_eigen(_start_ef + f, _start_es + s) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}

}

void SED_CM_FMM::fillSEDZDD(int &_f_sed, int &_s_sed, int &_start_tf, int &_end_tf, int &_start_ts, int &_end_ts, int &_start_ef, int &_end_ef, int &_start_es, int &_end_es)
{
	//#pragma omp parallel
	{
		int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
		VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
		VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3];
		Tetrahedron tet_plu, tet_min;
		value_t a_fld, a_src;//场与源的公共面面积
		value_t epsp, epsm;//源处介电常数
		int sedid_f, bedge_f, sedid_s, bedge_s;

		//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
		int vp_f, vp_s, vm_f, vm_s, v_ct1, v_ct2, v_ct3;

		Complex coef(0.0f, 0.0f);
		//Qcout << '\n' << "Fill the matrix ZDD:";
		//填充ZDD矩阵
		//tool::BarAndPercent bar_perc1;
		coef = Z0 * Z0 / (J0* PI4);
		//#pragma omp parallel
		for (int f = _start_tf; f < _end_tf; f++)
		{
			//bar_perc1(f-_start_tf + 1, _end_tf-_start_tf);

			ct_ptr_->getCommonTriangle(f, vp_f, vm_f, t_fld_plu, t_fld_min, a_fld, sedid_f, bedge_f);
			ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
			v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);

			shiftcoordv3(sedid_f, v_fld_ct3);//平移公共面三角形坐标

			v_fld_plu = mesh_ptr_->getVertex(vp_f);//场点+
			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);

			shiftcoordv1(sedid_f, v_fld_plu);//平移正顶点坐标
			shiftcoordv4(sedid_f, v_fld_plu4);//平移正四面体四点坐标

			int x1 = sedid_f / 3;
			int y1 = sedid_f % 3;

			if (vm_f != -1)
			{
				v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-
				tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-
				tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);


				shiftcoordv1(sedid_f, v_fld_min);
				shiftcoordv4(sedid_f, v_fld_min4);

				if (bedge_f != -1)//多一次平移边界处四面体
				{
					shiftbedgecoordv1(bedge_f, v_fld_min);
					shiftbedgecoordv4(bedge_f, v_fld_min4);
				}
			}

			for (int s = _start_ts; s < _end_ts; s++)
			{
				Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
				ct_ptr_->getCommonTriangle(s, vp_s, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm, sedid_s, bedge_s);
				ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
				v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
				v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
				v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);

				shiftcoordv3(sedid_s, v_src_ct3);

				v_src_plu = mesh_ptr_->getVertex(vp_s);
				tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
				tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);

				shiftcoordv1(sedid_s, v_src_plu);
				shiftcoordv4(sedid_s, v_src_plu4);

				int x2 = sedid_s / 3;
				int y2 = sedid_s % 3;

				if (vm_s != -1)
				{
					v_src_min = mesh_ptr_->getVertex(vm_s);
					tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);
					tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);

					shiftcoordv1(sedid_s, v_src_min);
					shiftcoordv4(sedid_s, v_src_min4);

					if (bedge_s != -1)
					{
						shiftbedgecoordv1(bedge_s, v_src_min);
						shiftbedgecoordv4(bedge_s, v_src_min4);
					}
				}
				//field tetrahedron+ <-->source tetrahedron+ 
				if ((t_fld_plu == t_src_plu) && (sedid_f == sedid_s))
				{
					Zpp = DDZppSingular(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
				}
				else
				{
					Zpp = DDZppKernel(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
				}
				//field tetrahedron+ <-->source tetrahedron-
				if (t_fld_plu == t_src_min)
				{
					if ((centerv4(v_fld_plu4) - centerv4(v_src_min4)).Norm() < 1e-5)
					{
						Zpm = DDZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
					}
					else
					{
						Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
					}
				}
				else if (t_src_min != -1)
				{
					Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
				}
				//field tetrahedron- <-->source tetrahedron+
				if (t_fld_min == t_src_plu)
				{
					if ((centerv4(v_fld_min4) - centerv4(v_src_plu4)).Norm() < 1e-5)
					{
						Zmp = DDZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
					}
					else
					{
						Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
					}
				}
				else if (t_fld_min != -1)
				{
					Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
				}
				//field tetrahedron- <-->source tetrahedron-
				if ((t_fld_min != -1) && (t_src_min != -1))
				{
					if (t_fld_min == t_src_min)
					{
						if ((centerv4(v_fld_min4) - centerv4(v_src_min4)).Norm() < 1e-5)
						{
							Zmm = DDZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
						}
						else
						{
							Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
						}
					}
					else
						Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
				}


				Z_sed(_start_ef + f, _start_es + s) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
				Z_sed_eigen(_start_ef + f, _start_es + s) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
			}
		}
	}

}

void SED_CM_FMM::fillSEDZDM(int &_f_sed, int &_s_sed, int &_start_tf, int &_end_tf, int &_start_ts, int &_end_ts, int &_start_ef, int &_end_ef, int &_start_es, int &_end_es)
{
	//#pragma omp parallel
	{
		int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
		int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle 
		VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
		VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
		Tetrahedron tet_plu, tet_min;
		Triangle tri_plu, tri_min;
		value_t a_fld, a_src;//SWG场与源的公共面面积
		value_t l_fld, l_src;//RWG场与源的公共边长度
		value_t epsp, epsm;//源处介电常数
		int sedid_f, bedge_f, sedid_s, bedge_s;

		//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
		int vp, vm_f, vm_s, vm, v_ct1, v_ct2, v_ct3;
		int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
		Complex coef(0.0f, 0.0f);
		//Qcout << '\n' << "Fill the matrix ZMD:";
		//填充ZDM矩阵
		//tool::BarAndPercent bar_perc2;
		coef = J0 * k_*Z0 / PI4;
		//#pragma omp parallel for
		for (int f = _start_tf; f < _end_tf; f++)
		{
			//bar_perc2(f - _start_tf + 1, _end_tf - _start_tf);
			ct_ptr_->getCommonTriangle(f, vp, vm_f, t_fld_plu, t_fld_min, a_fld, sedid_f, bedge_f);
			ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
			v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);

			shiftcoordv3(sedid_f, v_fld_ct3);

			v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);

			shiftcoordv1(sedid_f, v_fld_plu);
			shiftcoordv4(sedid_f, v_fld_plu4);

			if (vm_f != -1)
			{
				v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-
				tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-
				tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);

				shiftcoordv1(sedid_f, v_fld_min);
				shiftcoordv4(sedid_f, v_fld_min4);

				if (bedge_f != -1)
				{
					shiftbedgecoordv1(bedge_f, v_fld_min);
					shiftbedgecoordv4(bedge_f, v_fld_min4);
				}
			}
			for (int s = _start_es; s < _end_es; s++)
			{
				ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src, sedid_s, bedge_s);

				v_src_plu = mesh_ptr_->getVertex(vp);
				v_src_min = mesh_ptr_->getVertex(vm);

				tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
				tri_min = mesh_ptr_->getTriangleRef(f_src_min);

				tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
				tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);
				tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
				tri_min.getVertex(tri_min1, tri_min2, tri_min3);

				shiftcoordv1(sedid_s, v_src_plu);
				shiftcoordv1(sedid_s, v_src_min);
				shiftcoordv3(sedid_s, v_src_plu3);
				shiftcoordv3(sedid_s, v_src_min3);

				if (bedge_s != -1)
				{
					shiftbedgecoordv1(bedge_s, v_src_min);
					shiftbedgecoordv3(bedge_s, v_src_min3);
				}

				//fill
				Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);

				if (vm_f != -1)
				{
					Zmp = DMZmpKernel(v_fld_min4, v_src_plu3, v_fld_min, v_src_plu, v_fld_ct3, vm_f);
					Zmm = DMZmmKernel(v_fld_min4, v_src_min3, v_fld_min, v_src_min, v_fld_ct3, vm_f);
					Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
					Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
				}
				else if ((centerv3(v_fld_ct3) - centerv3(v_src_plu3)).Norm()<1e-5)
				{
					Zpp = DMZppSingular(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
					Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
				}
				else if ((centerv3(v_fld_ct3) - centerv3(v_src_min3)).Norm()<1e-5)
				{
					Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
					Zpm = DMZpmSingular(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
				}
				else
				{
					Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
					Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
				}

				Z_sed(_start_ef + f, _end_ts + s) = coef * a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
				Z_sed_eigen(_start_ef + f, _end_ts + s) = coef * a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
			}

		}
	}

}

void SED_CM_FMM::fillSEDZMD(int &_f_sed, int &_s_sed, int &_start_tf, int &_end_tf, int &_start_ts, int &_end_ts, int &_start_ef, int &_end_ef, int &_start_es, int &_end_es)
{
	//#pragma omp parallel
	{
		int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
		int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle 
		VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3], v_fld_plu3[3], v_fld_min3[3];
		VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3], v_src_plu3[3], v_src_min3[3];
		Tetrahedron tet_plu, tet_min;
		Triangle tri_plu, tri_min;
		value_t a_fld, a_src;//SWG场与源的公共面面积
		value_t l_fld, l_src;//RWG场与源的公共边长度
		value_t epsp, epsm;//源处介电常数
		int sedid_f, bedge_f, sedid_s, bedge_s;

		//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
		int vp, vm_f, vm_s, vm, v_ct1, v_ct2, v_ct3;
		int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
		Complex coef(0.0f, 0.0f);
		//Qcout << '\n' << "Fill the matrix ZDM:";
		//填充ZDM矩阵
		//tool::BarAndPercent bar_perc3;
		//coef = Z0 * Z0*k_ / (J0*PI4);
		//coef = Z0 * omiga*k_ / (PI4);
		coef = Z0 * k_ / (J0*PI4);
		//#pragma omp parallel for
		for (int f = _start_ef; f < _end_ef; f++)
		{
			//bar_perc3(f - _start_ef + 1, _end_ef - _start_ef);
			ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld, sedid_f, bedge_f);

			v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
			v_fld_min = mesh_ptr_->getVertex(vm);//场点-

			tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
			tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

			tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
			tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
			tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
			tri_min.getVertex(tri_min1, tri_min2, tri_min3);

			shiftcoordv1(sedid_f, v_fld_plu);
			shiftcoordv1(sedid_f, v_fld_min);
			shiftcoordv3(sedid_f, v_fld_plu3);
			shiftcoordv3(sedid_f, v_fld_min3);

			if (bedge_f != -1)
			{
				shiftbedgecoordv1(bedge_f, v_fld_min);
				shiftbedgecoordv3(bedge_f, v_fld_min3);
			}

			for (int s = _start_ts; s < _end_ts; s++)
			{
				Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
				ct_ptr_->getCommonTriangle(s, vp, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm, sedid_s, bedge_s);
				ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
				v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
				v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
				v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);

				shiftcoordv3(sedid_s, v_src_ct3);

				v_src_plu = mesh_ptr_->getVertex(vp);
				tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
				tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);

				shiftcoordv1(sedid_s, v_src_plu);
				shiftcoordv4(sedid_s, v_src_plu4);

				if (vm_s != -1)
				{
					v_src_min = mesh_ptr_->getVertex(vm_s);
					tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);
					tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);

					shiftcoordv1(sedid_s, v_src_min);
					shiftcoordv4(sedid_s, v_src_min4);

					if (bedge_s != -1)
					{
						shiftbedgecoordv1(bedge_s, v_src_min);
						shiftbedgecoordv4(bedge_s, v_src_min4);
					}
				}
				//fill
				if ((centerv3(v_fld_plu3) - centerv3(v_src_ct3)).Norm()<1e-5)
				{
					Zpp = MDZppSingular(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
				}
				else
					Zpp = MDZppkernel(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);


				if ((centerv3(v_fld_min3) - centerv3(v_src_ct3)).Norm()<1e-5)
				{
					Zmp = MDZmpSingular(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
				}
				else
					Zmp = MDZmpkernel(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);

				if (vm_s != -1)
				{
					Zpm = MDZpmkernel(v_fld_plu3, v_src_min4, v_fld_plu, v_src_min, v_src_ct3, vm_s, epsp, epsm);
					Zmm = MDZmmkernel(v_fld_min3, v_src_min4, v_fld_min, v_src_min, v_src_ct3, vm_s, epsp, epsm);
				}

				Z_sed(f + _end_tf, s + _start_es) = coef * l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
				Z_sed_eigen(f + _end_tf, s + _start_es) = coef * l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
			}
		}
	}

}

void SED_CM_FMM::fillSEDZMM(int &_f_sed, int &_s_sed, int &_start_tf, int &_end_tf, int &_start_ts, int &_end_ts, int &_start_ef, int &_end_ef, int &_start_es, int &_end_es)
{
	//#pragma omp parallel
	{
		int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle
		VectorR3 v_fld_plu, v_fld_min, v_fld_plu3[3], v_fld_min3[3];
		VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
		VectorR3 v_fld_cen, v_src_cen;
		Triangle tri_plu, tri_min;
		value_t l_fld, l_src;//场与源的公共边长度
		value_t scale_f, scale_s;//场源处的缩放因子
		int sedid_f, bedge_f, sedid_s, bedge_s;

		//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
		int vp, vm;
		Complex coef = (J0 * k_ * Z0) / PI4;
		//Qcout << '\n' << "Fill the matrix ZMM:";
		//填充ZMM矩阵
		//tool::BarAndPercent bar_perc4;
		//coef = (J0*k_*Z0) / PI4;
		//#pragma omp parallel for
		for (int f = _start_ef; f < _end_ef; f++)
		{
			//bar_perc4(f - _start_ef + 1, _end_ef - _start_ef);
			ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld, sedid_f, bedge_f);

			v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
			v_fld_min = mesh_ptr_->getVertex(vm);//场点-

			tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
			tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

			tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
			tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);

			/*shiftcoordv1(sedid_f, v_fld_plu);
			shiftcoordv1(sedid_f, v_fld_min);
			shiftcoordv3(sedid_f, v_fld_plu3);
			shiftcoordv3(sedid_f, v_fld_min3);*/

			shiftcoordv1(sedid_f, v_fld_plu, scale_x);
			shiftcoordv1(sedid_f, v_fld_min, scale_x);
			shiftcoordv3(sedid_f, v_fld_plu3, scale_x);
			shiftcoordv3(sedid_f, v_fld_min3, scale_x);
			scale_f = pow(scale_x, _f_sed / 3);
			l_fld *= scale_f;

			if (bedge_f != -1)
			{
				shiftbedgecoordv1(bedge_f, v_fld_min);
				shiftbedgecoordv3(bedge_f, v_fld_min3);
			}
			v_fld_cen = (centerv3(v_fld_plu3) + centerv3(v_fld_min3)) / 2.0f;

			for (int s = _start_es; s < _end_es; s++)
			{
				Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
				ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src, sedid_s, bedge_s);

				v_src_plu = mesh_ptr_->getVertex(vp);
				v_src_min = mesh_ptr_->getVertex(vm);

				tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
				tri_min = mesh_ptr_->getTriangleRef(f_src_min);

				tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
				tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);

				/*shiftcoordv1(sedid_s, v_src_plu);
				shiftcoordv1(sedid_s, v_src_min);
				shiftcoordv3(sedid_s, v_src_plu3);
				shiftcoordv3(sedid_s, v_src_min3);*/

				shiftcoordv1(sedid_s, v_src_plu, scale_x);
				shiftcoordv1(sedid_s, v_src_min, scale_x);
				shiftcoordv3(sedid_s, v_src_plu3, scale_x);
				shiftcoordv3(sedid_s, v_src_min3, scale_x);
				scale_s = pow(scale_x, _s_sed / 3);
				l_src *= scale_s;

				if (bedge_s != -1)
				{
					shiftbedgecoordv1(bedge_s, v_src_min);
					shiftbedgecoordv3(bedge_s, v_src_min3);
				}
				v_src_cen = (centerv3(v_src_plu3) + centerv3(v_src_min3)) / 2.0f;

				//////////filling////////////////////////////////
				if ((v_fld_cen - v_src_cen).Norm() > threshold_edm)
				{
					Z_sed_eigen(f + _end_tf, s + _end_ts) = MMEDMKernel(v_fld_plu3, v_fld_min3, v_src_plu3, v_src_min3, l_fld, l_src);
				}
				else
				{
					//field face+ <-->  source face+
					if ((f_fld_plu == f_src_plu) && (sedid_f == sedid_s))
					{
						Zpp = MMZppSingular(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
					}
					else
					{
						Zpp = MMZppKernel(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
					}
					//field face+ <--> source face-
					if ((f_fld_plu == f_src_min))
					{
						if ((centerv3(v_fld_plu3) - centerv3(v_src_min3)).Norm() < 1e-5)
						{
							Zpm = MMZpmSingular(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
						}
						else
						{
							Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
						}
					}
					else
					{
						Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
					}
					//field face- <--> source face+
					if ((f_fld_min == f_src_plu))
					{
						if ((centerv3(v_fld_min3) - centerv3(v_src_plu3)).Norm() < 1e-5)
						{
							Zmp = MMZmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
						}
						else
						{
							Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
						}
					}
					else
					{
						Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
					}
					//field face- <--> source face-
					if ((f_fld_min == f_src_min))
					{
						if ((centerv3(v_fld_min3) - centerv3(v_src_min3)).Norm() < 1e-5)
						{
							Zmm = MMZmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
						}
						else
						{
							Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
						}
					}
					else
					{
						Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
					}

					Z_sed(f + _end_tf, s + _end_ts) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
					Z_sed_eigen(f + _end_tf, s + _end_ts) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
				}
			}
		}
	}

}

void SED_CM_FMM::fillCOUPZDD(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n)
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3];
	Tetrahedron tet_plu, tet_min;
	value_t a_fld, a_src;//场与源的公共面面积
	value_t epsp, epsm;//源处介电常数

	int sedid_f, sedid_s;
	int bedge_f, bedge_s;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;
	int arrID_f, arrID_s;
	arrID_f = _x_m * Array_y + _y_m;
	arrID_s = _x_n * Array_y + _y_n;

	ct_ptr_->getSEDCommonTriangleSize(_k_m, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(_k_m, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(_k_n, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(_k_n, start_es, end_es);

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp_f, vp_s, vm_f, vm_s, v_ct1, v_ct2, v_ct3;

	Complex coef(0.0f, 0.0f);

	//填充ZDD矩阵
	//tool::BarAndPercent bar_perc1;
	coef = Z0 * Z0 / (J0* PI4);
	//#pragma omp parallel for
	for (int f = start_tf; f < end_tf; f++)
	{
		//(f - start_tf + 1, end_tf - start_tf);

		ct_ptr_->getCommonTriangle(f, vp_f, vm_f, t_fld_plu, t_fld_min, a_fld, sedid_f, bedge_f);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);

		shiftcoordv3(_x_m, _y_m, v_fld_ct3);

		v_fld_plu = mesh_ptr_->getVertex(vp_f);//场点+
		tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
		tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv4(_x_m, _y_m, v_fld_plu4);

		if (vm_f != -1)
		{

			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);


			shiftcoordv1(_x_m, _y_m, v_fld_min);
			shiftcoordv4(_x_m, _y_m, v_fld_min4);

			if (bedge_f != -1)
			{
				shiftbedgecoordv1(bedge_f, v_fld_min);
				shiftbedgecoordv4(bedge_f, v_fld_min4);
			}
		}


		for (int s = start_ts; s < end_ts; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp_s, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm, sedid_s, bedge_s);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);

			shiftcoordv3(_x_n, _y_n, v_src_ct3);

			v_src_plu = mesh_ptr_->getVertex(vp_s);
			tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
			tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv4(_x_n, _y_n, v_src_plu4);

			if (vm_s != -1)
			{
				v_src_min = mesh_ptr_->getVertex(vm_s);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);

				shiftcoordv1(_x_n, _y_n, v_src_min);
				shiftcoordv4(_x_n, _y_n, v_src_min4);

				if (bedge_s != -1)
				{
					shiftbedgecoordv1(bedge_s, v_src_min);
					shiftbedgecoordv4(bedge_s, v_src_min4);
				}
			}

			//field tetrahedron+ <-->source tetrahedron+ 
			if ((t_fld_plu == t_src_plu) && (arrID_f == arrID_s))
			{
				Zpp = DDZppSingular(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			}
			else
			{
				Zpp = DDZppKernel(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			}
			//field tetrahedron+ <-->source tetrahedron-
			if (t_fld_plu == t_src_min)
			{
				if ((centerv4(v_fld_plu4) - centerv4(v_src_min4)).Norm() < 1e-5)
				{
					Zpm = DDZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
				}
				else
				{
					Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
				}
			}
			else if (t_src_min != -1)
			{
				Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			}
			//field tetrahedron- <-->source tetrahedron+
			if (t_fld_min == t_src_plu)
			{
				if ((centerv4(v_fld_min4) - centerv4(v_src_plu4)).Norm() < 1e-5)
				{
					Zmp = DDZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
				}
				else
				{
					Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
				}
			}
			else if (t_fld_min != -1)
			{
				Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			}
			//field tetrahedron- <-->source tetrahedron-
			if ((t_fld_min != -1) && (t_src_min != -1))
			{
				if (t_fld_min == t_src_min)
				{
					if ((centerv4(v_fld_min4) - centerv4(v_src_min4)).Norm() < 1e-5)
					{
						Zmm = DDZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
					}
					else
					{
						Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
					}
				}
				else
					Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			}

			Z_coup(f - start_tf, s - start_ts) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
			Z_coup_eigen(f - start_tf, s - start_ts) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void SED_CM_FMM::fillCOUPZDM(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n)
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle 
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	value_t a_fld, a_src;//SWG场与源的公共面面积
	value_t l_fld, l_src;//RWG场与源的公共边长度
	value_t epsp, epsm;//源处介电常数

	int sedid_f, sedid_s;
	int bedge_f, bedge_s;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;

	/*int arrID_f, arrID_s;
	arrID_f = _x_m * Array_y + _y_m;
	arrID_s = _x_n * Array_y + _y_n;*/

	ct_ptr_->getSEDCommonTriangleSize(_k_m, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(_k_m, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(_k_n, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(_k_n, start_es, end_es);

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm_f, vm_s, vm, v_ct1, v_ct2, v_ct3;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f, 0.0f);

	//填充ZDM矩阵
	//tool::BarAndPercent bar_perc2;
	coef = J0 * k_*Z0 / PI4;
	//#pragma omp parallel for
	for (int f = start_tf; f < end_tf; f++)
	{
		//bar_perc2(f - start_tf + 1, end_tf - start_tf);
		ct_ptr_->getCommonTriangle(f, vp, vm_f, t_fld_plu, t_fld_min, a_fld, sedid_f, bedge_f);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);

		shiftcoordv3(_x_m, _y_m, v_fld_ct3);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
		tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv4(_x_m, _y_m, v_fld_plu4);

		if (vm_f != -1)
		{
			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);

			shiftcoordv1(_x_m, _y_m, v_fld_min);
			shiftcoordv4(_x_m, _y_m, v_fld_min4);

			if (bedge_f != -1)
			{
				shiftbedgecoordv1(bedge_f, v_fld_min);
				shiftbedgecoordv4(bedge_f, v_fld_min4);
			}
		}

		for (int s = start_es; s < end_es; s++)
		{
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src, sedid_s, bedge_s);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);
			tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
			tri_min.getVertex(tri_min1, tri_min2, tri_min3);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv1(_x_n, _y_n, v_src_min);
			shiftcoordv3(_x_n, _y_n, v_src_plu3);
			shiftcoordv3(_x_n, _y_n, v_src_min3);

			if (bedge_s != -1)
			{
				shiftbedgecoordv1(bedge_s, v_src_min);
				shiftbedgecoordv3(bedge_s, v_src_min3);
			}

			//fill
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);

			if (vm_f != -1)
			{
				Zmp = DMZmpKernel(v_fld_min4, v_src_plu3, v_fld_min, v_src_plu, v_fld_ct3, vm_f);
				Zmm = DMZmmKernel(v_fld_min4, v_src_min3, v_fld_min, v_src_min, v_fld_ct3, vm_f);
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else if ((centerv3(v_fld_ct3) - centerv3(v_src_plu3)).Norm()<1e-5)
			{
				Zpp = DMZppSingular(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else if ((centerv3(v_fld_ct3) - centerv3(v_src_min3)).Norm()<1e-5)
			{
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmSingular(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else
			{
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}

			Z_coup(f - start_tf, end_ts - start_ts + s - start_es) = coef * a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
			Z_coup_eigen(f - start_tf, end_ts - start_ts + s - start_es) = coef * a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
		}

	}
}

void SED_CM_FMM::fillCOUPZMD(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n)
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle 
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3], v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3], v_src_plu3[3], v_src_min3[3];
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	value_t a_fld, a_src;//SWG场与源的公共面面积
	value_t l_fld, l_src;//RWG场与源的公共边长度
	value_t epsp, epsm;//源处介电常数

	int sedid_f, sedid_s;
	int bedge_f, bedge_s;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;

	ct_ptr_->getSEDCommonTriangleSize(_k_m, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(_k_m, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(_k_n, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(_k_n, start_es, end_es);

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm_f, vm_s, vm, v_ct1, v_ct2, v_ct3;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f, 0.0f);

	//填充ZDM矩阵
	//tool::BarAndPercent bar_perc3;
	coef = Z0 * Z0*k_ / (J0*PI4);
	//#pragma omp parallel for
	for (int f = start_ef; f < end_ef; f++)
	{
		//bar_perc3(f - start_ef + 1, end_ef - start_ef);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld, sedid_f, bedge_f);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
		tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
		tri_min.getVertex(tri_min1, tri_min2, tri_min3);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv1(_x_m, _y_m, v_fld_min);
		shiftcoordv3(_x_m, _y_m, v_fld_plu3);
		shiftcoordv3(_x_m, _y_m, v_fld_min3);

		if (bedge_f != -1)
		{
			shiftbedgecoordv1(bedge_f, v_fld_min);
			shiftbedgecoordv3(bedge_f, v_fld_min3);
		}
		for (int s = start_ts; s < end_ts; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm, sedid_s, bedge_s);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);

			shiftcoordv3(_x_n, _y_n, v_src_ct3);

			v_src_plu = mesh_ptr_->getVertex(vp);
			tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
			tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv4(_x_n, _y_n, v_src_plu4);
			if (vm_s != -1)
			{
				v_src_min = mesh_ptr_->getVertex(vm_s);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);

				shiftcoordv1(_x_n, _y_n, v_src_min);
				shiftcoordv4(_x_n, _y_n, v_src_min4);

				if (bedge_s != -1)
				{
					shiftbedgecoordv1(bedge_s, v_src_min);
					shiftbedgecoordv4(bedge_s, v_src_min4);
				}
			}
			//fill
			if ((centerv3(v_fld_plu3) - centerv3(v_src_ct3)).Norm()<1e-5)
			{
				Zpp = MDZppSingular(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
				Zpp = MDZppkernel(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);


			if ((centerv3(v_fld_min3) - centerv3(v_src_ct3)).Norm()<1e-5)
			{
				Zmp = MDZmpSingular(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
				Zmp = MDZmpkernel(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);

			if (vm_s != -1)
			{
				Zpm = MDZpmkernel(v_fld_plu3, v_src_min4, v_fld_plu, v_src_min, v_src_ct3, vm_s, epsp, epsm);
				Zmm = MDZmmkernel(v_fld_min3, v_src_min4, v_fld_min, v_src_min, v_src_ct3, vm_s, epsp, epsm);
			}

			Z_coup(f + end_tf - start_tf - start_ef, s - start_ts) = coef * l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
			Z_coup_eigen(f + end_tf - start_tf - start_ef, s - start_ts) = coef * l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void SED_CM_FMM::fillCOUPZMM(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n)
{
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	Triangle tri_plu, tri_min;
	value_t l_fld, l_src;//场与源的公共边长度
	value_t scale_f, scale_s;//场源处的缩放因子

	int sedid_f, sedid_s;
	int bedge_f, bedge_s;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;
	int arrID_f, arrID_s;
	arrID_f = _x_m * Array_y + _y_m;
	arrID_s = _x_n * Array_y + _y_n;

	ct_ptr_->getSEDCommonTriangleSize(_k_m, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(_k_m, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(_k_n, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(_k_n, start_es, end_es);


	Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm;
	Complex coef = (J0 * k_ * Z0) / PI4;

	//填充ZMM矩阵
	//tool::BarAndPercent bar_perc4;
	//coef = (J0*k_*Z0) / PI4;
	//#pragma omp parallel for
	for (int f = start_ef; f < end_ef; f++)
	{
		//bar_perc4(f - start_ef + 1, end_ef - start_ef);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld, sedid_f, bedge_f);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu, scale_x);
		shiftcoordv1(_x_m, _y_m, v_fld_min, scale_x);
		shiftcoordv3(_x_m, _y_m, v_fld_plu3, scale_x);
		shiftcoordv3(_x_m, _y_m, v_fld_min3, scale_x);
		scale_f = pow(scale_x, _x_m);
		l_fld *= scale_f;

		if (bedge_f != -1)
		{
			shiftbedgecoordv1(bedge_f, v_fld_min);
			shiftbedgecoordv3(bedge_f, v_fld_min3);
		}
		for (int s = start_es; s < end_es; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src, sedid_s, bedge_s);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);

			shiftcoordv1(_x_n, _y_n, v_src_plu, scale_x);
			shiftcoordv1(_x_n, _y_n, v_src_min, scale_x);
			shiftcoordv3(_x_n, _y_n, v_src_plu3, scale_x);
			shiftcoordv3(_x_n, _y_n, v_src_min3, scale_x);
			scale_s = pow(scale_x, _x_n);
			l_src *= scale_s;

			if (bedge_s != -1)
			{
				shiftbedgecoordv1(bedge_s, v_src_min);
				shiftbedgecoordv3(bedge_s, v_src_min3);
			}
			//field face+ <-->  source face+
			if ((f_fld_plu == f_src_plu) && (arrID_f == arrID_s))
			{
				Zpp = MMZppSingular(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			else
			{
				Zpp = MMZppKernel(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			//field face+ <--> source face-
			if (f_fld_plu == f_src_min)
			{
				if ((centerv3(v_fld_plu3) - centerv3(v_src_min3)).Norm() < 1e-5)
				{
					Zpm = MMZpmSingular(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
				}
				else
				{
					Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
				}
			}
			else
			{
				Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
			}
			//field face- <--> source face+
			if (f_fld_min == f_src_plu)
			{
				if ((centerv3(v_fld_min3) - centerv3(v_src_plu3)).Norm() < 1e-5)
				{
					Zmp = MMZmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
				}
				else
				{
					Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
				}
			}
			else
			{
				Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
			}
			//field face- <--> source face-
			if (f_fld_min == f_src_min)
			{
				if ((centerv3(v_fld_min3) - centerv3(v_src_min3)).Norm() < 1e-5)
				{
					Zmm = MMZmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
				}
				else
				{
					Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
				}
			}
			else
			{
				Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
			}
			Z_coup(f + end_tf - start_tf - start_ef, s + end_ts - start_ts - start_es) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
			Z_coup_eigen(f + end_tf - start_tf - start_ef, s + end_ts - start_ts - start_es) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void SED_CM_FMM::fillCOUPZDD_MP(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup)
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3];
	Tetrahedron tet_plu, tet_min;
	value_t a_fld, a_src;//场与源的公共面面积
	value_t epsp, epsm;//源处介电常数

	int sedid_f, sedid_s;
	int bedge_f, bedge_s;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;
	int arrID_f, arrID_s;
	arrID_f = _x_m * Array_y + _y_m;
	arrID_s = _x_n * Array_y + _y_n;

	ct_ptr_->getSEDCommonTriangleSize(_k_m, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(_k_m, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(_k_n, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(_k_n, start_es, end_es);

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp_f, vp_s, vm_f, vm_s, v_ct1, v_ct2, v_ct3;

	Complex coef(0.0f, 0.0f);
	KIM ik = RED_KernelIntegral(_x_m, _y_m, _x_n, _y_n);
	//填充ZDD矩阵
	//tool::BarAndPercent bar_perc1;
	//coef = Z0 * Z0 / (J0* PI4);
	//coef = Z0 * omiga / (PI4);
	coef = Z0 / (J0* PI4);
	for (int f = start_tf; f < end_tf; f++)
	{
		//(f - start_tf + 1, end_tf - start_tf);

		ct_ptr_->getCommonTriangle(f, vp_f, vm_f, t_fld_plu, t_fld_min, a_fld, sedid_f, bedge_f);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);

		shiftcoordv3(_x_m, _y_m, v_fld_ct3);

		v_fld_plu = mesh_ptr_->getVertex(vp_f);//场点+
		tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
		tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv4(_x_m, _y_m, v_fld_plu4);

		if (vm_f != -1)
		{

			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);


			shiftcoordv1(_x_m, _y_m, v_fld_min);
			shiftcoordv4(_x_m, _y_m, v_fld_min4);

			if (bedge_f != -1)
			{
				shiftbedgecoordv1(bedge_f, v_fld_min);
				shiftbedgecoordv4(bedge_f, v_fld_min4);
			}
		}


		for (int s = start_ts; s < end_ts; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp_s, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm, sedid_s, bedge_s);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);

			shiftcoordv3(_x_n, _y_n, v_src_ct3);

			v_src_plu = mesh_ptr_->getVertex(vp_s);
			tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
			tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv4(_x_n, _y_n, v_src_plu4);

			if (vm_s != -1)
			{
				v_src_min = mesh_ptr_->getVertex(vm_s);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);

				shiftcoordv1(_x_n, _y_n, v_src_min);
				shiftcoordv4(_x_n, _y_n, v_src_min4);

				if (bedge_s != -1)
				{
					shiftbedgecoordv1(bedge_s, v_src_min);
					shiftbedgecoordv4(bedge_s, v_src_min4);
				}
			}
			if (bedge_f == -1 && bedge_s == -1)
			{
				if (arrID_f == arrID_s)
				{
					if (t_fld_plu == t_src_plu)
						Zpp = DDZppSingular_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
					else
						Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_plu == t_src_min)
						Zpm = DDZpmSingular_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);
					else if (t_src_min != -1)
						Zpm = DDZpmKernel_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_min == t_src_plu)
						Zmp = DDZmpSingular_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);
					else if (t_fld_min != -1)
						Zmp = DDZmpKernel_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_min != -1 && t_src_min != -1)
					{
						if (t_fld_min == t_src_min)
							Zmm = DDZmmSingular_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
						else
							Zmm = DDZmmKernel_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
					}

				}
				else
				{
					Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);

					if (t_src_min != -1)
						Zpm = DDZpmKernel_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);

					if (t_fld_min != -1)
						Zmp = DDZmpKernel_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);

					if (t_fld_min != -1 && t_src_min != -1)
						Zmm = DDZmmKernel_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
				}
			}
			else
			{
				//field tetrahedron+ <-->source tetrahedron+ 
				if ((t_fld_plu == t_src_plu) && (arrID_f == arrID_s))
				{
					Zpp = DDZppSingular_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
				}
				else
				{
					Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
				}
				//field tetrahedron+ <-->source tetrahedron-
				if (t_fld_plu == t_src_min)
				{
					if ((centerv4(v_fld_plu4) - centerv4(v_src_min4)).Norm() < 1e-5)
					{
						Zpm = DDZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
					}
					else
					{
						Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
					}
				}
				else if (t_src_min != -1)
				{
					Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
				}
				//field tetrahedron- <-->source tetrahedron+
				if (t_fld_min == t_src_plu)
				{
					if ((centerv4(v_fld_min4) - centerv4(v_src_plu4)).Norm() < 1e-5)
					{
						Zmp = DDZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
					}
					else
					{
						Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
					}
				}
				else if (t_fld_min != -1)
				{
					Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
				}
				//field tetrahedron- <-->source tetrahedron-
				if ((t_fld_min != -1) && (t_src_min != -1))
				{
					if (t_fld_min == t_src_min)
					{
						if ((centerv4(v_fld_min4) - centerv4(v_src_min4)).Norm() < 1e-5)
						{
							Zmm = DDZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
						}
						else
						{
							Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
						}
					}
					else
						Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
				}
			}
			////field tetrahedron+ <-->source tetrahedron+ 
			//if ((t_fld_plu == t_src_plu) && (arrID_f == arrID_s))
			//{
			//	Zpp = DDZppSingular(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			//}
			//else
			//{
			//	Zpp = DDZppKernel(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			//}
			////field tetrahedron+ <-->source tetrahedron-
			//if (t_fld_plu == t_src_min)
			//{
			//	if ((centerv4(v_fld_plu4) - centerv4(v_src_min4)).Norm() < 1e-5)
			//	{
			//		Zpm = DDZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			//	}
			//	else
			//	{
			//		Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			//	}
			//}
			//else if (t_src_min != -1)
			//{
			//	Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			//}
			////field tetrahedron- <-->source tetrahedron+
			//if (t_fld_min == t_src_plu)
			//{
			//	if ((centerv4(v_fld_min4) - centerv4(v_src_plu4)).Norm() < 1e-5)
			//	{
			//		Zmp = DDZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			//	}
			//	else
			//	{
			//		Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			//	}
			//}
			//else if (t_fld_min != -1)
			//{
			//	Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			//}
			////field tetrahedron- <-->source tetrahedron-
			//if ((t_fld_min != -1) && (t_src_min != -1))
			//{
			//	if (t_fld_min == t_src_min)
			//	{
			//		if ((centerv4(v_fld_min4) - centerv4(v_src_min4)).Norm() < 1e-5)
			//		{
			//			Zmm = DDZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			//		}
			//		else
			//		{
			//			Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			//		}
			//	}
			//	else
			//		Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			//}

			Z_coup(f - start_tf, s - start_ts) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f - start_tf, s - start_ts) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void SED_CM_FMM::fillCOUPZDM_MP(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup)
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle 
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	value_t a_fld, a_src;//SWG场与源的公共面面积
	value_t l_fld, l_src;//RWG场与源的公共边长度
	value_t epsp, epsm;//源处介电常数

	int sedid_f, sedid_s;
	int bedge_f, bedge_s;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;

	/*int arrID_f, arrID_s;
	arrID_f = _x_m * Array_y + _y_m;
	arrID_s = _x_n * Array_y + _y_n;*/

	ct_ptr_->getSEDCommonTriangleSize(_k_m, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(_k_m, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(_k_n, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(_k_n, start_es, end_es);

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm_f, vm_s, vm, v_ct1, v_ct2, v_ct3;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f, 0.0f);

	//填充ZDM矩阵
	//tool::BarAndPercent bar_perc2;
	coef = J0 * k_*Z0 / PI4;
	//#pragma omp parallel for
	for (int f = start_tf; f < end_tf; f++)
	{
		//bar_perc2(f - start_tf + 1, end_tf - start_tf);
		ct_ptr_->getCommonTriangle(f, vp, vm_f, t_fld_plu, t_fld_min, a_fld, sedid_f, bedge_f);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);

		shiftcoordv3(_x_m, _y_m, v_fld_ct3);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
		tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv4(_x_m, _y_m, v_fld_plu4);

		if (vm_f != -1)
		{
			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);

			shiftcoordv1(_x_m, _y_m, v_fld_min);
			shiftcoordv4(_x_m, _y_m, v_fld_min4);

			if (bedge_f != -1)
			{
				shiftbedgecoordv1(bedge_f, v_fld_min);
				shiftbedgecoordv4(bedge_f, v_fld_min4);
			}
		}

		for (int s = start_es; s < end_es; s++)
		{
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src, sedid_s, bedge_s);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);
			tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
			tri_min.getVertex(tri_min1, tri_min2, tri_min3);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv1(_x_n, _y_n, v_src_min);
			shiftcoordv3(_x_n, _y_n, v_src_plu3);
			shiftcoordv3(_x_n, _y_n, v_src_min3);

			if (bedge_s != -1)
			{
				shiftbedgecoordv1(bedge_s, v_src_min);
				shiftbedgecoordv3(bedge_s, v_src_min3);
			}

			//fill
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);

			if (vm_f != -1)
			{
				Zmp = DMZmpKernel(v_fld_min4, v_src_plu3, v_fld_min, v_src_plu, v_fld_ct3, vm_f);
				Zmm = DMZmmKernel(v_fld_min4, v_src_min3, v_fld_min, v_src_min, v_fld_ct3, vm_f);
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else if ((centerv3(v_fld_ct3) - centerv3(v_src_plu3)).Norm()<1e-5)
			{
				Zpp = DMZppSingular(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else if ((centerv3(v_fld_ct3) - centerv3(v_src_min3)).Norm()<1e-5)
			{
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmSingular(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else
			{
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}

			Z_coup(f - start_tf, end_ts - start_ts + s - start_es) = coef * a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f - start_tf, end_ts - start_ts + s - start_es) = coef * a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
		}

	}
}

void SED_CM_FMM::fillCOUPZMD_MP(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup)
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle 
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3], v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3], v_src_plu3[3], v_src_min3[3];
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	value_t a_fld, a_src;//SWG场与源的公共面面积
	value_t l_fld, l_src;//RWG场与源的公共边长度
	value_t epsp, epsm;//源处介电常数

	int sedid_f, sedid_s;
	int bedge_f, bedge_s;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;

	ct_ptr_->getSEDCommonTriangleSize(_k_m, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(_k_m, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(_k_n, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(_k_n, start_es, end_es);

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm_f, vm_s, vm, v_ct1, v_ct2, v_ct3;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f, 0.0f);

	//填充ZDM矩阵
	//tool::BarAndPercent bar_perc3;
	//coef = Z0 * Z0*k_ / (J0*PI4);
	//coef = Z0 * omiga*k_ / (PI4);
	coef = Z0 * k_ / (J0*PI4);
	for (int f = start_ef; f < end_ef; f++)
	{
		//bar_perc3(f - start_ef + 1, end_ef - start_ef);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld, sedid_f, bedge_f);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
		tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
		tri_min.getVertex(tri_min1, tri_min2, tri_min3);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv1(_x_m, _y_m, v_fld_min);
		shiftcoordv3(_x_m, _y_m, v_fld_plu3);
		shiftcoordv3(_x_m, _y_m, v_fld_min3);

		if (bedge_f != -1)
		{
			shiftbedgecoordv1(bedge_f, v_fld_min);
			shiftbedgecoordv3(bedge_f, v_fld_min3);
		}
		for (int s = start_ts; s < end_ts; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm, sedid_s, bedge_s);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);

			shiftcoordv3(_x_n, _y_n, v_src_ct3);

			v_src_plu = mesh_ptr_->getVertex(vp);
			tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
			tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv4(_x_n, _y_n, v_src_plu4);
			if (vm_s != -1)
			{
				v_src_min = mesh_ptr_->getVertex(vm_s);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);

				shiftcoordv1(_x_n, _y_n, v_src_min);
				shiftcoordv4(_x_n, _y_n, v_src_min4);

				if (bedge_s != -1)
				{
					shiftbedgecoordv1(bedge_s, v_src_min);
					shiftbedgecoordv4(bedge_s, v_src_min4);
				}
			}
			//fill
			if ((centerv3(v_fld_plu3) - centerv3(v_src_ct3)).Norm()<1e-5)
			{
				Zpp = MDZppSingular(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
				Zpp = MDZppkernel(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);


			if ((centerv3(v_fld_min3) - centerv3(v_src_ct3)).Norm()<1e-5)
			{
				Zmp = MDZmpSingular(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
				Zmp = MDZmpkernel(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);

			if (vm_s != -1)
			{
				Zpm = MDZpmkernel(v_fld_plu3, v_src_min4, v_fld_plu, v_src_min, v_src_ct3, vm_s, epsp, epsm);
				Zmm = MDZmmkernel(v_fld_min3, v_src_min4, v_fld_min, v_src_min, v_src_ct3, vm_s, epsp, epsm);
			}

			Z_coup(f + end_tf - start_tf - start_ef, s - start_ts) = coef * l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f + end_tf - start_tf - start_ef, s - start_ts) = coef * l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void SED_CM_FMM::fillCOUPZMM_MP(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup)
{
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	VectorR3 v_fld_cen, v_src_cen;
	Triangle tri_plu, tri_min;
	value_t l_fld, l_src;//场与源的公共边长度

	int sedid_f, sedid_s;
	int bedge_f, bedge_s;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;
	int arrID_f, arrID_s;
	arrID_f = _x_m * Array_y + _y_m;
	arrID_s = _x_n * Array_y + _y_n;

	ct_ptr_->getSEDCommonTriangleSize(_k_m, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(_k_m, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(_k_n, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(_k_n, start_es, end_es);


	Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm;
	Complex coef = (J0 * k_ * Z0) / PI4;

	//填充ZMM矩阵
	//tool::BarAndPercent bar_perc4;
	//coef = (J0*k_*Z0) / PI4;
	//#pragma omp parallel for
	for (int f = start_ef; f < end_ef; f++)
	{
		//bar_perc4(f - start_ef + 1, end_ef - start_ef);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld, sedid_f, bedge_f);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu, scale_x);
		shiftcoordv1(_x_m, _y_m, v_fld_min, scale_x);
		shiftcoordv3(_x_m, _y_m, v_fld_plu3, scale_x);
		shiftcoordv3(_x_m, _y_m, v_fld_min3, scale_x);
		l_fld *= pow(scale_x, _x_m);
		//Qcout << _x_m << _y_m << std::endl;
		//system("pause");
		if (bedge_f != -1)
		{
			shiftbedgecoordv1(bedge_f, v_fld_min);
			shiftbedgecoordv3(bedge_f, v_fld_min3);
		}
		v_fld_cen = (centerv3(v_fld_plu3) + centerv3(v_fld_min3)) / 2.0f;

		for (int s = start_es; s < end_es; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src, sedid_s, bedge_s);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);

			shiftcoordv1(_x_n, _y_n, v_src_plu, scale_x);
			shiftcoordv1(_x_n, _y_n, v_src_min, scale_x);
			shiftcoordv3(_x_n, _y_n, v_src_plu3, scale_x);
			shiftcoordv3(_x_n, _y_n, v_src_min3, scale_x);
			l_src *= pow(scale_x, _x_n);
			//Qcout << _x_n << _y_n << std::endl;
			//system("pause");
			if (bedge_s != -1)
			{
				shiftbedgecoordv1(bedge_s, v_src_min);
				shiftbedgecoordv3(bedge_s, v_src_min3);
			}
			v_src_cen = (centerv3(v_src_plu3) + centerv3(v_src_min3)) / 2.0f;

			//////////////////////filling/////////////////////////////////////////////////
			if ((v_fld_cen - v_src_cen).Norm() > threshold_edm)
			{
				Z_coup(f + end_tf - start_tf - start_ef, s + end_ts - start_ts - start_es) = MMEDMKernel(v_fld_plu3, v_fld_min3, v_src_plu3, v_src_min3, l_fld, l_src);
			}
			else
			{
				//field face+ <-->  source face+
				if ((f_fld_plu == f_src_plu) && (arrID_f == arrID_s))
				{
					Zpp = MMZppSingular(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
				}
				else
				{
					Zpp = MMZppKernel(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
				}
				//field face+ <--> source face-
				if (f_fld_plu == f_src_min)
				{
					if ((centerv3(v_fld_plu3) - centerv3(v_src_min3)).Norm() < 1e-5)
					{
						Zpm = MMZpmSingular(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
					}
					else
					{
						Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
					}
				}
				else
				{
					Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
				}
				//field face- <--> source face+
				if (f_fld_min == f_src_plu)
				{
					if ((centerv3(v_fld_min3) - centerv3(v_src_plu3)).Norm() < 1e-5)
					{
						Zmp = MMZmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
					}
					else
					{
						Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
					}
				}
				else
				{
					Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
				}
				//field face- <--> source face-
				if (f_fld_min == f_src_min)
				{
					if ((centerv3(v_fld_min3) - centerv3(v_src_min3)).Norm() < 1e-5)
					{
						Zmm = MMZmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
					}
					else
					{
						Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
					}
				}
				else
				{
					Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
				}
				Z_coup(f + end_tf - start_tf - start_ef, s + end_ts - start_ts - start_es) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
			}
		}
	}
}

void SED_CM_FMM::fillCOUPZDD_MP_C(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup)
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3];
	Tetrahedron tet_plu, tet_min;
	value_t a_fld, a_src;//场与源的公共面面积
	value_t epsp, epsm;//源处介电常数

	int sedid_f, sedid_s;
	int bedge_f, bedge_s;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;
	int arrID_f, arrID_s;
	arrID_f = _x_m * Array_y + _y_m;
	arrID_s = _x_n * Array_y + _y_n;

	ct_ptr_->getSEDCommonTriangleSize(_k_m, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(_k_m, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(_k_n, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(_k_n, start_es, end_es);

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp_f, vp_s, vm_f, vm_s, v_ct1, v_ct2, v_ct3;

	Complex coef(0.0f, 0.0f);
	KIM ik = RED_KernelIntegral(_x_m, _y_m, _x_n, _y_n);
	//填充ZDD矩阵
	//tool::BarAndPercent bar_perc1;
	//coef = Z0 * Z0 / (J0* PI4);
	//coef = Z0 * omiga / (PI4);
	coef = Z0 / (J0* PI4);
	for (int f = start_tf; f < start_tf + com_ct; f++)
	{
		//(f - start_tf + 1, end_tf - start_tf);

		ct_ptr_->getCommonTriangle(f, vp_f, vm_f, t_fld_plu, t_fld_min, a_fld, sedid_f, bedge_f);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);

		shiftcoordv3(_x_m, _y_m, v_fld_ct3);

		v_fld_plu = mesh_ptr_->getVertex(vp_f);//场点+
		tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
		tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv4(_x_m, _y_m, v_fld_plu4);

		if (vm_f != -1)
		{

			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);


			shiftcoordv1(_x_m, _y_m, v_fld_min);
			shiftcoordv4(_x_m, _y_m, v_fld_min4);

			if (bedge_f != -1)
			{
				shiftbedgecoordv1(bedge_f, v_fld_min);
				shiftbedgecoordv4(bedge_f, v_fld_min4);
			}
		}


		for (int s = start_ts; s < start_ts + com_ct; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp_s, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm, sedid_s, bedge_s);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);

			shiftcoordv3(_x_n, _y_n, v_src_ct3);

			v_src_plu = mesh_ptr_->getVertex(vp_s);
			tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
			tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv4(_x_n, _y_n, v_src_plu4);

			if (vm_s != -1)
			{
				v_src_min = mesh_ptr_->getVertex(vm_s);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);

				shiftcoordv1(_x_n, _y_n, v_src_min);
				shiftcoordv4(_x_n, _y_n, v_src_min4);

				if (bedge_s != -1)
				{
					shiftbedgecoordv1(bedge_s, v_src_min);
					shiftbedgecoordv4(bedge_s, v_src_min4);
				}
			}
			if (bedge_f == -1 && bedge_s == -1)
			{
				if (arrID_f == arrID_s)
				{
					if (t_fld_plu == t_src_plu)
						Zpp = DDZppSingular_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
					else
						Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_plu == t_src_min)
						Zpm = DDZpmSingular_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);
					else if (t_src_min != -1)
						Zpm = DDZpmKernel_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_min == t_src_plu)
						Zmp = DDZmpSingular_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);
					else if (t_fld_min != -1)
						Zmp = DDZmpKernel_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_min != -1 && t_src_min != -1)
					{
						if (t_fld_min == t_src_min)
							Zmm = DDZmmSingular_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
						else
							Zmm = DDZmmKernel_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
					}

				}
				else
				{
					Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);

					if (t_src_min != -1)
						Zpm = DDZpmKernel_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);

					if (t_fld_min != -1)
						Zmp = DDZmpKernel_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);

					if (t_fld_min != -1 && t_src_min != -1)
						Zmm = DDZmmKernel_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
				}
			}
			else
			{
				//field tetrahedron+ <-->source tetrahedron+ 
				if ((t_fld_plu == t_src_plu) && (arrID_f == arrID_s))
				{
					Zpp = DDZppSingular_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
				}
				else
				{
					Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
				}
				//field tetrahedron+ <-->source tetrahedron-
				if (t_fld_plu == t_src_min)
				{
					if ((centerv4(v_fld_plu4) - centerv4(v_src_min4)).Norm() < 1e-5)
					{
						Zpm = DDZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
					}
					else
					{
						Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
					}
				}
				else if (t_src_min != -1)
				{
					Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
				}
				//field tetrahedron- <-->source tetrahedron+
				if (t_fld_min == t_src_plu)
				{
					if ((centerv4(v_fld_min4) - centerv4(v_src_plu4)).Norm() < 1e-5)
					{
						Zmp = DDZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
					}
					else
					{
						Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
					}
				}
				else if (t_fld_min != -1)
				{
					Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
				}
				//field tetrahedron- <-->source tetrahedron-
				if ((t_fld_min != -1) && (t_src_min != -1))
				{
					if (t_fld_min == t_src_min)
					{
						if ((centerv4(v_fld_min4) - centerv4(v_src_min4)).Norm() < 1e-5)
						{
							Zmm = DDZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
						}
						else
						{
							Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
						}
					}
					else
						Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
				}
			}
			////field tetrahedron+ <-->source tetrahedron+ 
			//if ((t_fld_plu == t_src_plu) && (arrID_f == arrID_s))
			//{
			//	Zpp = DDZppSingular(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			//}
			//else
			//{
			//	Zpp = DDZppKernel(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			//}
			////field tetrahedron+ <-->source tetrahedron-
			//if (t_fld_plu == t_src_min)
			//{
			//	if ((centerv4(v_fld_plu4) - centerv4(v_src_min4)).Norm() < 1e-5)
			//	{
			//		Zpm = DDZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			//	}
			//	else
			//	{
			//		Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			//	}
			//}
			//else if (t_src_min != -1)
			//{
			//	Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			//}
			////field tetrahedron- <-->source tetrahedron+
			//if (t_fld_min == t_src_plu)
			//{
			//	if ((centerv4(v_fld_min4) - centerv4(v_src_plu4)).Norm() < 1e-5)
			//	{
			//		Zmp = DDZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			//	}
			//	else
			//	{
			//		Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			//	}
			//}
			//else if (t_fld_min != -1)
			//{
			//	Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			//}
			////field tetrahedron- <-->source tetrahedron-
			//if ((t_fld_min != -1) && (t_src_min != -1))
			//{
			//	if (t_fld_min == t_src_min)
			//	{
			//		if ((centerv4(v_fld_min4) - centerv4(v_src_min4)).Norm() < 1e-5)
			//		{
			//			Zmm = DDZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			//		}
			//		else
			//		{
			//			Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			//		}
			//	}
			//	else
			//		Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			//}

			Z_coup(f - start_tf, s - start_ts) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f - start_tf, s - start_ts) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void SED_CM_FMM::fillCOUPZDM_MP_C(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup)
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle 
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	value_t a_fld, a_src;//SWG场与源的公共面面积
	value_t l_fld, l_src;//RWG场与源的公共边长度
	value_t epsp, epsm;//源处介电常数

	int sedid_f, sedid_s;
	int bedge_f, bedge_s;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;

	/*int arrID_f, arrID_s;
	arrID_f = _x_m * Array_y + _y_m;
	arrID_s = _x_n * Array_y + _y_n;*/

	ct_ptr_->getSEDCommonTriangleSize(_k_m, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(_k_m, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(_k_n, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(_k_n, start_es, end_es);

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm_f, vm_s, vm, v_ct1, v_ct2, v_ct3;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f, 0.0f);

	//填充ZDM矩阵
	//tool::BarAndPercent bar_perc2;
	coef = J0 * k_*Z0 / PI4;
	//#pragma omp parallel for
	for (int f = start_tf; f < start_tf + com_ct; f++)
	{
		//bar_perc2(f - start_tf + 1, end_tf - start_tf);
		ct_ptr_->getCommonTriangle(f, vp, vm_f, t_fld_plu, t_fld_min, a_fld, sedid_f, bedge_f);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);

		shiftcoordv3(_x_m, _y_m, v_fld_ct3);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
		tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv4(_x_m, _y_m, v_fld_plu4);

		if (vm_f != -1)
		{
			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);

			shiftcoordv1(_x_m, _y_m, v_fld_min);
			shiftcoordv4(_x_m, _y_m, v_fld_min4);

			if (bedge_f != -1)
			{
				shiftbedgecoordv1(bedge_f, v_fld_min);
				shiftbedgecoordv4(bedge_f, v_fld_min4);
			}
		}

		for (int s = start_es; s < start_es + com_ce; s++)
		{
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src, sedid_s, bedge_s);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);
			tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
			tri_min.getVertex(tri_min1, tri_min2, tri_min3);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv1(_x_n, _y_n, v_src_min);
			shiftcoordv3(_x_n, _y_n, v_src_plu3);
			shiftcoordv3(_x_n, _y_n, v_src_min3);

			if (bedge_s != -1)
			{
				shiftbedgecoordv1(bedge_s, v_src_min);
				shiftbedgecoordv3(bedge_s, v_src_min3);
			}

			//fill
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);

			if (vm_f != -1)
			{
				Zmp = DMZmpKernel(v_fld_min4, v_src_plu3, v_fld_min, v_src_plu, v_fld_ct3, vm_f);
				Zmm = DMZmmKernel(v_fld_min4, v_src_min3, v_fld_min, v_src_min, v_fld_ct3, vm_f);
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else if ((centerv3(v_fld_ct3) - centerv3(v_src_plu3)).Norm()<1e-5)
			{
				Zpp = DMZppSingular(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else if ((centerv3(v_fld_ct3) - centerv3(v_src_min3)).Norm()<1e-5)
			{
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmSingular(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else
			{
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}

			Z_coup(f - start_tf, com_ct + s - start_es) = coef * a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f - start_tf, end_ts - start_ts + s - start_es) = coef * a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
		}

	}
}

void SED_CM_FMM::fillCOUPZMD_MP_C(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup)
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle 
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3], v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3], v_src_plu3[3], v_src_min3[3];
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	value_t a_fld, a_src;//SWG场与源的公共面面积
	value_t l_fld, l_src;//RWG场与源的公共边长度
	value_t epsp, epsm;//源处介电常数

	int sedid_f, sedid_s;
	int bedge_f, bedge_s;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;

	ct_ptr_->getSEDCommonTriangleSize(_k_m, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(_k_m, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(_k_n, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(_k_n, start_es, end_es);

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm_f, vm_s, vm, v_ct1, v_ct2, v_ct3;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f, 0.0f);

	//填充ZDM矩阵
	//tool::BarAndPercent bar_perc3;
	//coef = Z0 * Z0*k_ / (J0*PI4);
	//coef = Z0 * omiga*k_ / (PI4);
	coef = Z0 * k_ / (J0*PI4);
	for (int f = start_ef; f < start_ef + com_ce; f++)
	{
		//bar_perc3(f - start_ef + 1, end_ef - start_ef);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld, sedid_f, bedge_f);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
		tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
		tri_min.getVertex(tri_min1, tri_min2, tri_min3);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv1(_x_m, _y_m, v_fld_min);
		shiftcoordv3(_x_m, _y_m, v_fld_plu3);
		shiftcoordv3(_x_m, _y_m, v_fld_min3);

		if (bedge_f != -1)
		{
			shiftbedgecoordv1(bedge_f, v_fld_min);
			shiftbedgecoordv3(bedge_f, v_fld_min3);
		}
		for (int s = start_ts; s < start_ts + com_ct; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm, sedid_s, bedge_s);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);

			shiftcoordv3(_x_n, _y_n, v_src_ct3);

			v_src_plu = mesh_ptr_->getVertex(vp);
			tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
			tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv4(_x_n, _y_n, v_src_plu4);
			if (vm_s != -1)
			{
				v_src_min = mesh_ptr_->getVertex(vm_s);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);

				shiftcoordv1(_x_n, _y_n, v_src_min);
				shiftcoordv4(_x_n, _y_n, v_src_min4);

				if (bedge_s != -1)
				{
					shiftbedgecoordv1(bedge_s, v_src_min);
					shiftbedgecoordv4(bedge_s, v_src_min4);
				}
			}
			//fill
			if ((centerv3(v_fld_plu3) - centerv3(v_src_ct3)).Norm()<1e-5)
			{
				Zpp = MDZppSingular(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
				Zpp = MDZppkernel(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);


			if ((centerv3(v_fld_min3) - centerv3(v_src_ct3)).Norm()<1e-5)
			{
				Zmp = MDZmpSingular(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
				Zmp = MDZmpkernel(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);

			if (vm_s != -1)
			{
				Zpm = MDZpmkernel(v_fld_plu3, v_src_min4, v_fld_plu, v_src_min, v_src_ct3, vm_s, epsp, epsm);
				Zmm = MDZmmkernel(v_fld_min3, v_src_min4, v_fld_min, v_src_min, v_src_ct3, vm_s, epsp, epsm);
			}

			Z_coup(com_ct + f - start_ef, s - start_ts) = coef * l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f + end_tf - start_tf - start_ef, s - start_ts) = coef * l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void SED_CM_FMM::fillCOUPZMM_MP_C(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup)
{
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	Triangle tri_plu, tri_min;
	value_t l_fld, l_src;//场与源的公共边长度

	int sedid_f, sedid_s;
	int bedge_f, bedge_s;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;
	int arrID_f, arrID_s;
	arrID_f = _x_m * Array_y + _y_m;
	arrID_s = _x_n * Array_y + _y_n;

	ct_ptr_->getSEDCommonTriangleSize(_k_m, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(_k_m, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(_k_n, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(_k_n, start_es, end_es);


	Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm;
	Complex coef = (J0 * k_ * Z0) / PI4;

	//填充ZMM矩阵
	//tool::BarAndPercent bar_perc4;
	//coef = (J0*k_*Z0) / PI4;
	//#pragma omp parallel for
	for (int f = start_ef; f < start_ef + com_ce; f++)
	{
		//bar_perc4(f - start_ef + 1, end_ef - start_ef);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld, sedid_f, bedge_f);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv1(_x_m, _y_m, v_fld_min);
		shiftcoordv3(_x_m, _y_m, v_fld_plu3);
		shiftcoordv3(_x_m, _y_m, v_fld_min3);

		if (bedge_f != -1)
		{
			shiftbedgecoordv1(bedge_f, v_fld_min);
			shiftbedgecoordv3(bedge_f, v_fld_min3);
		}
		for (int s = start_es; s < start_es + com_ce; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src, sedid_s, bedge_s);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv1(_x_n, _y_n, v_src_min);
			shiftcoordv3(_x_n, _y_n, v_src_plu3);
			shiftcoordv3(_x_n, _y_n, v_src_min3);

			if (bedge_s != -1)
			{
				shiftbedgecoordv1(bedge_s, v_src_min);
				shiftbedgecoordv3(bedge_s, v_src_min3);
			}
			//field face+ <-->  source face+
			if ((f_fld_plu == f_src_plu) && (arrID_f == arrID_s))
			{
				Zpp = MMZppSingular(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			else
			{
				Zpp = MMZppKernel(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			//field face+ <--> source face-
			if (f_fld_plu == f_src_min)
			{
				if ((centerv3(v_fld_plu3) - centerv3(v_src_min3)).Norm() < 1e-5)
				{
					Zpm = MMZpmSingular(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
				}
				else
				{
					Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
				}
			}
			else
			{
				Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
			}
			//field face- <--> source face+
			if (f_fld_min == f_src_plu)
			{
				if ((centerv3(v_fld_min3) - centerv3(v_src_plu3)).Norm() < 1e-5)
				{
					Zmp = MMZmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
				}
				else
				{
					Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
				}
			}
			else
			{
				Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
			}
			//field face- <--> source face-
			if (f_fld_min == f_src_min)
			{
				if ((centerv3(v_fld_min3) - centerv3(v_src_min3)).Norm() < 1e-5)
				{
					Zmm = MMZmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
				}
				else
				{
					Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
				}
			}
			else
			{
				Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
			}
			Z_coup(f + com_ct - start_ef, s + com_ct - start_es) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f + end_tf - start_tf - start_ef, s + end_ts - start_ts - start_es) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void SED_CM_FMM::fillCOUPZDD_MP_B(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup)
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3];
	Tetrahedron tet_plu, tet_min;
	value_t a_fld, a_src;//场与源的公共面面积
	value_t epsp, epsm;//源处介电常数

	int sedid_f, sedid_s;
	int bedge_f, bedge_s;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;
	int arrID_f, arrID_s;
	arrID_f = _x_m * Array_y + _y_m;
	arrID_s = _x_n * Array_y + _y_n;

	ct_ptr_->getSEDCommonTriangleSize(_k_m, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(_k_m, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(_k_n, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(_k_n, start_es, end_es);

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp_f, vp_s, vm_f, vm_s, v_ct1, v_ct2, v_ct3;

	Complex coef(0.0f, 0.0f);
	KIM ik = RED_KernelIntegral(_x_m, _y_m, _x_n, _y_n);
	//填充ZDD矩阵
	//tool::BarAndPercent bar_perc1;
	//coef = Z0 * Z0 / (J0* PI4);
	//coef = Z0 * omiga / (PI4);
	coef = Z0 / (J0* PI4);
	for (int f = start_tf; f < end_tf; f++)
	{
		//(f - start_tf + 1, end_tf - start_tf);

		ct_ptr_->getCommonTriangle(f, vp_f, vm_f, t_fld_plu, t_fld_min, a_fld, sedid_f, bedge_f);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);

		shiftcoordv3(_x_m, _y_m, v_fld_ct3);

		v_fld_plu = mesh_ptr_->getVertex(vp_f);//场点+
		tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
		tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv4(_x_m, _y_m, v_fld_plu4);

		if (vm_f != -1)
		{

			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);


			shiftcoordv1(_x_m, _y_m, v_fld_min);
			shiftcoordv4(_x_m, _y_m, v_fld_min4);

			if (bedge_f != -1)
			{
				shiftbedgecoordv1(bedge_f, v_fld_min);
				shiftbedgecoordv4(bedge_f, v_fld_min4);
			}
		}


		for (int s = start_ts + com_ct; s < end_ts; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp_s, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm, sedid_s, bedge_s);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);

			shiftcoordv3(_x_n, _y_n, v_src_ct3);

			v_src_plu = mesh_ptr_->getVertex(vp_s);
			tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
			tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv4(_x_n, _y_n, v_src_plu4);

			if (vm_s != -1)
			{
				v_src_min = mesh_ptr_->getVertex(vm_s);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);

				shiftcoordv1(_x_n, _y_n, v_src_min);
				shiftcoordv4(_x_n, _y_n, v_src_min4);

				if (bedge_s != -1)
				{
					shiftbedgecoordv1(bedge_s, v_src_min);
					shiftbedgecoordv4(bedge_s, v_src_min4);
				}
			}
			if (bedge_f == -1 && bedge_s == -1)
			{
				if (arrID_f == arrID_s)
				{
					if (t_fld_plu == t_src_plu)
						Zpp = DDZppSingular_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
					else
						Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_plu == t_src_min)
						Zpm = DDZpmSingular_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);
					else if (t_src_min != -1)
						Zpm = DDZpmKernel_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_min == t_src_plu)
						Zmp = DDZmpSingular_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);
					else if (t_fld_min != -1)
						Zmp = DDZmpKernel_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_min != -1 && t_src_min != -1)
					{
						if (t_fld_min == t_src_min)
							Zmm = DDZmmSingular_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
						else
							Zmm = DDZmmKernel_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
					}

				}
				else
				{
					Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);

					if (t_src_min != -1)
						Zpm = DDZpmKernel_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);

					if (t_fld_min != -1)
						Zmp = DDZmpKernel_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);

					if (t_fld_min != -1 && t_src_min != -1)
						Zmm = DDZmmKernel_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
				}
			}
			else
			{
				//field tetrahedron+ <-->source tetrahedron+ 
				if ((t_fld_plu == t_src_plu) && (arrID_f == arrID_s))
				{
					Zpp = DDZppSingular_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
				}
				else
				{
					Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
				}
				//field tetrahedron+ <-->source tetrahedron-
				if (t_fld_plu == t_src_min)
				{
					if ((centerv4(v_fld_plu4) - centerv4(v_src_min4)).Norm() < 1e-5)
					{
						Zpm = DDZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
					}
					else
					{
						Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
					}
				}
				else if (t_src_min != -1)
				{
					Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
				}
				//field tetrahedron- <-->source tetrahedron+
				if (t_fld_min == t_src_plu)
				{
					if ((centerv4(v_fld_min4) - centerv4(v_src_plu4)).Norm() < 1e-5)
					{
						Zmp = DDZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
					}
					else
					{
						Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
					}
				}
				else if (t_fld_min != -1)
				{
					Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
				}
				//field tetrahedron- <-->source tetrahedron-
				if ((t_fld_min != -1) && (t_src_min != -1))
				{
					if (t_fld_min == t_src_min)
					{
						if ((centerv4(v_fld_min4) - centerv4(v_src_min4)).Norm() < 1e-5)
						{
							Zmm = DDZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
						}
						else
						{
							Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
						}
					}
					else
						Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
				}
			}
			////field tetrahedron+ <-->source tetrahedron+ 
			//if ((t_fld_plu == t_src_plu) && (arrID_f == arrID_s))
			//{
			//	Zpp = DDZppSingular(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			//}
			//else
			//{
			//	Zpp = DDZppKernel(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			//}
			////field tetrahedron+ <-->source tetrahedron-
			//if (t_fld_plu == t_src_min)
			//{
			//	if ((centerv4(v_fld_plu4) - centerv4(v_src_min4)).Norm() < 1e-5)
			//	{
			//		Zpm = DDZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			//	}
			//	else
			//	{
			//		Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			//	}
			//}
			//else if (t_src_min != -1)
			//{
			//	Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			//}
			////field tetrahedron- <-->source tetrahedron+
			//if (t_fld_min == t_src_plu)
			//{
			//	if ((centerv4(v_fld_min4) - centerv4(v_src_plu4)).Norm() < 1e-5)
			//	{
			//		Zmp = DDZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			//	}
			//	else
			//	{
			//		Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			//	}
			//}
			//else if (t_fld_min != -1)
			//{
			//	Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			//}
			////field tetrahedron- <-->source tetrahedron-
			//if ((t_fld_min != -1) && (t_src_min != -1))
			//{
			//	if (t_fld_min == t_src_min)
			//	{
			//		if ((centerv4(v_fld_min4) - centerv4(v_src_min4)).Norm() < 1e-5)
			//		{
			//			Zmm = DDZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			//		}
			//		else
			//		{
			//			Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			//		}
			//	}
			//	else
			//		Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			//}

			Z_coup(f - start_tf, s - start_ts) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f - start_tf, s - start_ts) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
	for (int f = start_tf + com_ct; f < end_tf; f++)
	{
		//(f - start_tf + 1, end_tf - start_tf);

		ct_ptr_->getCommonTriangle(f, vp_f, vm_f, t_fld_plu, t_fld_min, a_fld, sedid_f, bedge_f);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);

		shiftcoordv3(_x_m, _y_m, v_fld_ct3);

		v_fld_plu = mesh_ptr_->getVertex(vp_f);//场点+
		tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
		tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv4(_x_m, _y_m, v_fld_plu4);

		if (vm_f != -1)
		{

			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);


			shiftcoordv1(_x_m, _y_m, v_fld_min);
			shiftcoordv4(_x_m, _y_m, v_fld_min4);

			if (bedge_f != -1)
			{
				shiftbedgecoordv1(bedge_f, v_fld_min);
				shiftbedgecoordv4(bedge_f, v_fld_min4);
			}
		}


		for (int s = start_ts; s < end_ts; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp_s, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm, sedid_s, bedge_s);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);

			shiftcoordv3(_x_n, _y_n, v_src_ct3);

			v_src_plu = mesh_ptr_->getVertex(vp_s);
			tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
			tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv4(_x_n, _y_n, v_src_plu4);

			if (vm_s != -1)
			{
				v_src_min = mesh_ptr_->getVertex(vm_s);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);

				shiftcoordv1(_x_n, _y_n, v_src_min);
				shiftcoordv4(_x_n, _y_n, v_src_min4);

				if (bedge_s != -1)
				{
					shiftbedgecoordv1(bedge_s, v_src_min);
					shiftbedgecoordv4(bedge_s, v_src_min4);
				}
			}
			if (bedge_f == -1 && bedge_s == -1)
			{
				if (arrID_f == arrID_s)
				{
					if (t_fld_plu == t_src_plu)
						Zpp = DDZppSingular_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
					else
						Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_plu == t_src_min)
						Zpm = DDZpmSingular_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);
					else if (t_src_min != -1)
						Zpm = DDZpmKernel_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_min == t_src_plu)
						Zmp = DDZmpSingular_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);
					else if (t_fld_min != -1)
						Zmp = DDZmpKernel_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_min != -1 && t_src_min != -1)
					{
						if (t_fld_min == t_src_min)
							Zmm = DDZmmSingular_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
						else
							Zmm = DDZmmKernel_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
					}

				}
				else
				{
					Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);

					if (t_src_min != -1)
						Zpm = DDZpmKernel_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);

					if (t_fld_min != -1)
						Zmp = DDZmpKernel_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);

					if (t_fld_min != -1 && t_src_min != -1)
						Zmm = DDZmmKernel_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
				}
			}
			else
			{
				//field tetrahedron+ <-->source tetrahedron+ 
				if ((t_fld_plu == t_src_plu) && (arrID_f == arrID_s))
				{
					Zpp = DDZppSingular_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
				}
				else
				{
					Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
				}
				//field tetrahedron+ <-->source tetrahedron-
				if (t_fld_plu == t_src_min)
				{
					if ((centerv4(v_fld_plu4) - centerv4(v_src_min4)).Norm() < 1e-5)
					{
						Zpm = DDZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
					}
					else
					{
						Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
					}
				}
				else if (t_src_min != -1)
				{
					Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
				}
				//field tetrahedron- <-->source tetrahedron+
				if (t_fld_min == t_src_plu)
				{
					if ((centerv4(v_fld_min4) - centerv4(v_src_plu4)).Norm() < 1e-5)
					{
						Zmp = DDZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
					}
					else
					{
						Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
					}
				}
				else if (t_fld_min != -1)
				{
					Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
				}
				//field tetrahedron- <-->source tetrahedron-
				if ((t_fld_min != -1) && (t_src_min != -1))
				{
					if (t_fld_min == t_src_min)
					{
						if ((centerv4(v_fld_min4) - centerv4(v_src_min4)).Norm() < 1e-5)
						{
							Zmm = DDZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
						}
						else
						{
							Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
						}
					}
					else
						Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
				}
			}
			////field tetrahedron+ <-->source tetrahedron+ 
			//if ((t_fld_plu == t_src_plu) && (arrID_f == arrID_s))
			//{
			//	Zpp = DDZppSingular(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			//}
			//else
			//{
			//	Zpp = DDZppKernel(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			//}
			////field tetrahedron+ <-->source tetrahedron-
			//if (t_fld_plu == t_src_min)
			//{
			//	if ((centerv4(v_fld_plu4) - centerv4(v_src_min4)).Norm() < 1e-5)
			//	{
			//		Zpm = DDZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			//	}
			//	else
			//	{
			//		Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			//	}
			//}
			//else if (t_src_min != -1)
			//{
			//	Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			//}
			////field tetrahedron- <-->source tetrahedron+
			//if (t_fld_min == t_src_plu)
			//{
			//	if ((centerv4(v_fld_min4) - centerv4(v_src_plu4)).Norm() < 1e-5)
			//	{
			//		Zmp = DDZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			//	}
			//	else
			//	{
			//		Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			//	}
			//}
			//else if (t_fld_min != -1)
			//{
			//	Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			//}
			////field tetrahedron- <-->source tetrahedron-
			//if ((t_fld_min != -1) && (t_src_min != -1))
			//{
			//	if (t_fld_min == t_src_min)
			//	{
			//		if ((centerv4(v_fld_min4) - centerv4(v_src_min4)).Norm() < 1e-5)
			//		{
			//			Zmm = DDZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			//		}
			//		else
			//		{
			//			Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			//		}
			//	}
			//	else
			//		Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			//}

			Z_coup(f - start_tf, s - start_ts) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f - start_tf, s - start_ts) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
	for (int f = start_tf + com_ct; f < end_tf; f++)
	{
		//(f - start_tf + 1, end_tf - start_tf);

		ct_ptr_->getCommonTriangle(f, vp_f, vm_f, t_fld_plu, t_fld_min, a_fld, sedid_f, bedge_f);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);

		shiftcoordv3(_x_m, _y_m, v_fld_ct3);

		v_fld_plu = mesh_ptr_->getVertex(vp_f);//场点+
		tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
		tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv4(_x_m, _y_m, v_fld_plu4);

		if (vm_f != -1)
		{

			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);


			shiftcoordv1(_x_m, _y_m, v_fld_min);
			shiftcoordv4(_x_m, _y_m, v_fld_min4);

			if (bedge_f != -1)
			{
				shiftbedgecoordv1(bedge_f, v_fld_min);
				shiftbedgecoordv4(bedge_f, v_fld_min4);
			}
		}


		for (int s = start_ts + com_ct; s < end_ts; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp_s, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm, sedid_s, bedge_s);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);

			shiftcoordv3(_x_n, _y_n, v_src_ct3);

			v_src_plu = mesh_ptr_->getVertex(vp_s);
			tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
			tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv4(_x_n, _y_n, v_src_plu4);

			if (vm_s != -1)
			{
				v_src_min = mesh_ptr_->getVertex(vm_s);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);

				shiftcoordv1(_x_n, _y_n, v_src_min);
				shiftcoordv4(_x_n, _y_n, v_src_min4);

				if (bedge_s != -1)
				{
					shiftbedgecoordv1(bedge_s, v_src_min);
					shiftbedgecoordv4(bedge_s, v_src_min4);
				}
			}
			if (bedge_f == -1 && bedge_s == -1)
			{
				if (arrID_f == arrID_s)
				{
					if (t_fld_plu == t_src_plu)
						Zpp = DDZppSingular_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
					else
						Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_plu == t_src_min)
						Zpm = DDZpmSingular_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);
					else if (t_src_min != -1)
						Zpm = DDZpmKernel_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_min == t_src_plu)
						Zmp = DDZmpSingular_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);
					else if (t_fld_min != -1)
						Zmp = DDZmpKernel_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);
					///////////////////////////////////
					if (t_fld_min != -1 && t_src_min != -1)
					{
						if (t_fld_min == t_src_min)
							Zmm = DDZmmSingular_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
						else
							Zmm = DDZmmKernel_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
					}

				}
				else
				{
					Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);

					if (t_src_min != -1)
						Zpm = DDZpmKernel_fast(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm, ik);

					if (t_fld_min != -1)
						Zmp = DDZmpKernel_fast(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm, ik);

					if (t_fld_min != -1 && t_src_min != -1)
						Zmm = DDZmmKernel_fast(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm, ik);
				}
			}
			else
			{
				//field tetrahedron+ <-->source tetrahedron+ 
				if ((t_fld_plu == t_src_plu) && (arrID_f == arrID_s))
				{
					Zpp = DDZppSingular_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
				}
				else
				{
					Zpp = DDZppKernel_fast(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm, ik);
				}
				//field tetrahedron+ <-->source tetrahedron-
				if (t_fld_plu == t_src_min)
				{
					if ((centerv4(v_fld_plu4) - centerv4(v_src_min4)).Norm() < 1e-5)
					{
						Zpm = DDZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
					}
					else
					{
						Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
					}
				}
				else if (t_src_min != -1)
				{
					Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
				}
				//field tetrahedron- <-->source tetrahedron+
				if (t_fld_min == t_src_plu)
				{
					if ((centerv4(v_fld_min4) - centerv4(v_src_plu4)).Norm() < 1e-5)
					{
						Zmp = DDZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
					}
					else
					{
						Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
					}
				}
				else if (t_fld_min != -1)
				{
					Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
				}
				//field tetrahedron- <-->source tetrahedron-
				if ((t_fld_min != -1) && (t_src_min != -1))
				{
					if (t_fld_min == t_src_min)
					{
						if ((centerv4(v_fld_min4) - centerv4(v_src_min4)).Norm() < 1e-5)
						{
							Zmm = DDZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
						}
						else
						{
							Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
						}
					}
					else
						Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
				}
			}
			////field tetrahedron+ <-->source tetrahedron+ 
			//if ((t_fld_plu == t_src_plu) && (arrID_f == arrID_s))
			//{
			//	Zpp = DDZppSingular(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			//}
			//else
			//{
			//	Zpp = DDZppKernel(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			//}
			////field tetrahedron+ <-->source tetrahedron-
			//if (t_fld_plu == t_src_min)
			//{
			//	if ((centerv4(v_fld_plu4) - centerv4(v_src_min4)).Norm() < 1e-5)
			//	{
			//		Zpm = DDZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			//	}
			//	else
			//	{
			//		Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			//	}
			//}
			//else if (t_src_min != -1)
			//{
			//	Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			//}
			////field tetrahedron- <-->source tetrahedron+
			//if (t_fld_min == t_src_plu)
			//{
			//	if ((centerv4(v_fld_min4) - centerv4(v_src_plu4)).Norm() < 1e-5)
			//	{
			//		Zmp = DDZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			//	}
			//	else
			//	{
			//		Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			//	}
			//}
			//else if (t_fld_min != -1)
			//{
			//	Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			//}
			////field tetrahedron- <-->source tetrahedron-
			//if ((t_fld_min != -1) && (t_src_min != -1))
			//{
			//	if (t_fld_min == t_src_min)
			//	{
			//		if ((centerv4(v_fld_min4) - centerv4(v_src_min4)).Norm() < 1e-5)
			//		{
			//			Zmm = DDZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			//		}
			//		else
			//		{
			//			Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			//		}
			//	}
			//	else
			//		Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			//}

			Z_coup(f - start_tf, s - start_ts) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f - start_tf, s - start_ts) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void SED_CM_FMM::fillCOUPZDM_MP_B(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup)
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle 
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	value_t a_fld, a_src;//SWG场与源的公共面面积
	value_t l_fld, l_src;//RWG场与源的公共边长度
	value_t epsp, epsm;//源处介电常数

	int sedid_f, sedid_s;
	int bedge_f, bedge_s;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;

	/*int arrID_f, arrID_s;
	arrID_f = _x_m * Array_y + _y_m;
	arrID_s = _x_n * Array_y + _y_n;*/

	ct_ptr_->getSEDCommonTriangleSize(_k_m, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(_k_m, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(_k_n, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(_k_n, start_es, end_es);

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm_f, vm_s, vm, v_ct1, v_ct2, v_ct3;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f, 0.0f);

	//填充ZDM矩阵
	//tool::BarAndPercent bar_perc2;
	coef = J0 * k_*Z0 / PI4;
	//#pragma omp parallel for
	for (int f = start_tf; f < end_tf; f++)
	{
		//bar_perc2(f - start_tf + 1, end_tf - start_tf);
		ct_ptr_->getCommonTriangle(f, vp, vm_f, t_fld_plu, t_fld_min, a_fld, sedid_f, bedge_f);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);

		shiftcoordv3(_x_m, _y_m, v_fld_ct3);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
		tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv4(_x_m, _y_m, v_fld_plu4);

		if (vm_f != -1)
		{
			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);

			shiftcoordv1(_x_m, _y_m, v_fld_min);
			shiftcoordv4(_x_m, _y_m, v_fld_min4);

			if (bedge_f != -1)
			{
				shiftbedgecoordv1(bedge_f, v_fld_min);
				shiftbedgecoordv4(bedge_f, v_fld_min4);
			}
		}

		for (int s = start_es + com_ce; s < end_es; s++)
		{
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src, sedid_s, bedge_s);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);
			tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
			tri_min.getVertex(tri_min1, tri_min2, tri_min3);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv1(_x_n, _y_n, v_src_min);
			shiftcoordv3(_x_n, _y_n, v_src_plu3);
			shiftcoordv3(_x_n, _y_n, v_src_min3);

			if (bedge_s != -1)
			{
				shiftbedgecoordv1(bedge_s, v_src_min);
				shiftbedgecoordv3(bedge_s, v_src_min3);
			}

			//fill
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);

			if (vm_f != -1)
			{
				Zmp = DMZmpKernel(v_fld_min4, v_src_plu3, v_fld_min, v_src_plu, v_fld_ct3, vm_f);
				Zmm = DMZmmKernel(v_fld_min4, v_src_min3, v_fld_min, v_src_min, v_fld_ct3, vm_f);
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else if ((centerv3(v_fld_ct3) - centerv3(v_src_plu3)).Norm()<1e-5)
			{
				Zpp = DMZppSingular(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else if ((centerv3(v_fld_ct3) - centerv3(v_src_min3)).Norm()<1e-5)
			{
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmSingular(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else
			{
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}

			Z_coup(f - start_tf, end_ts - start_ts + s - start_es) = coef * a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f - start_tf, end_ts - start_ts + s - start_es) = coef * a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
		}

	}
	for (int f = start_tf + com_ct; f < end_tf; f++)
	{
		//bar_perc2(f - start_tf + 1, end_tf - start_tf);
		ct_ptr_->getCommonTriangle(f, vp, vm_f, t_fld_plu, t_fld_min, a_fld, sedid_f, bedge_f);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);

		shiftcoordv3(_x_m, _y_m, v_fld_ct3);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
		tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv4(_x_m, _y_m, v_fld_plu4);

		if (vm_f != -1)
		{
			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);

			shiftcoordv1(_x_m, _y_m, v_fld_min);
			shiftcoordv4(_x_m, _y_m, v_fld_min4);

			if (bedge_f != -1)
			{
				shiftbedgecoordv1(bedge_f, v_fld_min);
				shiftbedgecoordv4(bedge_f, v_fld_min4);
			}
		}

		for (int s = start_es; s < end_es; s++)
		{
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src, sedid_s, bedge_s);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);
			tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
			tri_min.getVertex(tri_min1, tri_min2, tri_min3);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv1(_x_n, _y_n, v_src_min);
			shiftcoordv3(_x_n, _y_n, v_src_plu3);
			shiftcoordv3(_x_n, _y_n, v_src_min3);

			if (bedge_s != -1)
			{
				shiftbedgecoordv1(bedge_s, v_src_min);
				shiftbedgecoordv3(bedge_s, v_src_min3);
			}

			//fill
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);

			if (vm_f != -1)
			{
				Zmp = DMZmpKernel(v_fld_min4, v_src_plu3, v_fld_min, v_src_plu, v_fld_ct3, vm_f);
				Zmm = DMZmmKernel(v_fld_min4, v_src_min3, v_fld_min, v_src_min, v_fld_ct3, vm_f);
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else if ((centerv3(v_fld_ct3) - centerv3(v_src_plu3)).Norm()<1e-5)
			{
				Zpp = DMZppSingular(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else if ((centerv3(v_fld_ct3) - centerv3(v_src_min3)).Norm()<1e-5)
			{
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmSingular(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else
			{
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}

			Z_coup(f - start_tf, end_ts - start_ts + s - start_es) = coef * a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f - start_tf, end_ts - start_ts + s - start_es) = coef * a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
		}

	}
	for (int f = start_tf + com_ct; f < end_tf; f++)
	{
		//bar_perc2(f - start_tf + 1, end_tf - start_tf);
		ct_ptr_->getCommonTriangle(f, vp, vm_f, t_fld_plu, t_fld_min, a_fld, sedid_f, bedge_f);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);

		shiftcoordv3(_x_m, _y_m, v_fld_ct3);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
		tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv4(_x_m, _y_m, v_fld_plu4);

		if (vm_f != -1)
		{
			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);

			shiftcoordv1(_x_m, _y_m, v_fld_min);
			shiftcoordv4(_x_m, _y_m, v_fld_min4);

			if (bedge_f != -1)
			{
				shiftbedgecoordv1(bedge_f, v_fld_min);
				shiftbedgecoordv4(bedge_f, v_fld_min4);
			}
		}

		for (int s = start_es + com_ce; s < end_es; s++)
		{
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src, sedid_s, bedge_s);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);
			tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
			tri_min.getVertex(tri_min1, tri_min2, tri_min3);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv1(_x_n, _y_n, v_src_min);
			shiftcoordv3(_x_n, _y_n, v_src_plu3);
			shiftcoordv3(_x_n, _y_n, v_src_min3);

			if (bedge_s != -1)
			{
				shiftbedgecoordv1(bedge_s, v_src_min);
				shiftbedgecoordv3(bedge_s, v_src_min3);
			}

			//fill
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);

			if (vm_f != -1)
			{
				Zmp = DMZmpKernel(v_fld_min4, v_src_plu3, v_fld_min, v_src_plu, v_fld_ct3, vm_f);
				Zmm = DMZmmKernel(v_fld_min4, v_src_min3, v_fld_min, v_src_min, v_fld_ct3, vm_f);
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else if ((centerv3(v_fld_ct3) - centerv3(v_src_plu3)).Norm()<1e-5)
			{
				Zpp = DMZppSingular(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else if ((centerv3(v_fld_ct3) - centerv3(v_src_min3)).Norm()<1e-5)
			{
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmSingular(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else
			{
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}

			Z_coup(f - start_tf, end_ts - start_ts + s - start_es) = coef * a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f - start_tf, end_ts - start_ts + s - start_es) = coef * a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
		}

	}
}

void SED_CM_FMM::fillCOUPZMD_MP_B(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup)
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle 
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3], v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3], v_src_plu3[3], v_src_min3[3];
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	value_t a_fld, a_src;//SWG场与源的公共面面积
	value_t l_fld, l_src;//RWG场与源的公共边长度
	value_t epsp, epsm;//源处介电常数

	int sedid_f, sedid_s;
	int bedge_f, bedge_s;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;

	ct_ptr_->getSEDCommonTriangleSize(_k_m, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(_k_m, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(_k_n, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(_k_n, start_es, end_es);

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm_f, vm_s, vm, v_ct1, v_ct2, v_ct3;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f, 0.0f);

	//填充ZDM矩阵
	//tool::BarAndPercent bar_perc3;
	//coef = Z0 * Z0*k_ / (J0*PI4);
	//coef = Z0 * omiga*k_ / (PI4);
	coef = Z0 * k_ / (J0*PI4);
	for (int f = start_ef; f < end_ef; f++)
	{
		//bar_perc3(f - start_ef + 1, end_ef - start_ef);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld, sedid_f, bedge_f);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
		tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
		tri_min.getVertex(tri_min1, tri_min2, tri_min3);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv1(_x_m, _y_m, v_fld_min);
		shiftcoordv3(_x_m, _y_m, v_fld_plu3);
		shiftcoordv3(_x_m, _y_m, v_fld_min3);

		if (bedge_f != -1)
		{
			shiftbedgecoordv1(bedge_f, v_fld_min);
			shiftbedgecoordv3(bedge_f, v_fld_min3);
		}
		for (int s = start_ts + com_ct; s < end_ts; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm, sedid_s, bedge_s);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);

			shiftcoordv3(_x_n, _y_n, v_src_ct3);

			v_src_plu = mesh_ptr_->getVertex(vp);
			tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
			tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv4(_x_n, _y_n, v_src_plu4);
			if (vm_s != -1)
			{
				v_src_min = mesh_ptr_->getVertex(vm_s);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);

				shiftcoordv1(_x_n, _y_n, v_src_min);
				shiftcoordv4(_x_n, _y_n, v_src_min4);

				if (bedge_s != -1)
				{
					shiftbedgecoordv1(bedge_s, v_src_min);
					shiftbedgecoordv4(bedge_s, v_src_min4);
				}
			}
			//fill
			if ((centerv3(v_fld_plu3) - centerv3(v_src_ct3)).Norm()<1e-5)
			{
				Zpp = MDZppSingular(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
				Zpp = MDZppkernel(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);


			if ((centerv3(v_fld_min3) - centerv3(v_src_ct3)).Norm()<1e-5)
			{
				Zmp = MDZmpSingular(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
				Zmp = MDZmpkernel(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);

			if (vm_s != -1)
			{
				Zpm = MDZpmkernel(v_fld_plu3, v_src_min4, v_fld_plu, v_src_min, v_src_ct3, vm_s, epsp, epsm);
				Zmm = MDZmmkernel(v_fld_min3, v_src_min4, v_fld_min, v_src_min, v_src_ct3, vm_s, epsp, epsm);
			}

			Z_coup(f + end_tf - start_tf - start_ef, s - start_ts) = coef * l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f + end_tf - start_tf - start_ef, s - start_ts) = coef * l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
		}
	}
	for (int f = start_ef + com_ce; f < end_ef; f++)
	{
		//bar_perc3(f - start_ef + 1, end_ef - start_ef);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld, sedid_f, bedge_f);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
		tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
		tri_min.getVertex(tri_min1, tri_min2, tri_min3);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv1(_x_m, _y_m, v_fld_min);
		shiftcoordv3(_x_m, _y_m, v_fld_plu3);
		shiftcoordv3(_x_m, _y_m, v_fld_min3);

		if (bedge_f != -1)
		{
			shiftbedgecoordv1(bedge_f, v_fld_min);
			shiftbedgecoordv3(bedge_f, v_fld_min3);
		}
		for (int s = start_ts; s < end_ts; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm, sedid_s, bedge_s);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);

			shiftcoordv3(_x_n, _y_n, v_src_ct3);

			v_src_plu = mesh_ptr_->getVertex(vp);
			tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
			tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv4(_x_n, _y_n, v_src_plu4);
			if (vm_s != -1)
			{
				v_src_min = mesh_ptr_->getVertex(vm_s);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);

				shiftcoordv1(_x_n, _y_n, v_src_min);
				shiftcoordv4(_x_n, _y_n, v_src_min4);

				if (bedge_s != -1)
				{
					shiftbedgecoordv1(bedge_s, v_src_min);
					shiftbedgecoordv4(bedge_s, v_src_min4);
				}
			}
			//fill
			if ((centerv3(v_fld_plu3) - centerv3(v_src_ct3)).Norm()<1e-5)
			{
				Zpp = MDZppSingular(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
				Zpp = MDZppkernel(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);


			if ((centerv3(v_fld_min3) - centerv3(v_src_ct3)).Norm()<1e-5)
			{
				Zmp = MDZmpSingular(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
				Zmp = MDZmpkernel(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);

			if (vm_s != -1)
			{
				Zpm = MDZpmkernel(v_fld_plu3, v_src_min4, v_fld_plu, v_src_min, v_src_ct3, vm_s, epsp, epsm);
				Zmm = MDZmmkernel(v_fld_min3, v_src_min4, v_fld_min, v_src_min, v_src_ct3, vm_s, epsp, epsm);
			}

			Z_coup(f + end_tf - start_tf - start_ef, s - start_ts) = coef * l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f + end_tf - start_tf - start_ef, s - start_ts) = coef * l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
		}
	}
	for (int f = start_ef + com_ce; f < end_ef; f++)
	{
		//bar_perc3(f - start_ef + 1, end_ef - start_ef);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld, sedid_f, bedge_f);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
		tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
		tri_min.getVertex(tri_min1, tri_min2, tri_min3);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv1(_x_m, _y_m, v_fld_min);
		shiftcoordv3(_x_m, _y_m, v_fld_plu3);
		shiftcoordv3(_x_m, _y_m, v_fld_min3);

		if (bedge_f != -1)
		{
			shiftbedgecoordv1(bedge_f, v_fld_min);
			shiftbedgecoordv3(bedge_f, v_fld_min3);
		}
		for (int s = start_ts + com_ct; s < end_ts; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm, sedid_s, bedge_s);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);

			shiftcoordv3(_x_n, _y_n, v_src_ct3);

			v_src_plu = mesh_ptr_->getVertex(vp);
			tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
			tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv4(_x_n, _y_n, v_src_plu4);
			if (vm_s != -1)
			{
				v_src_min = mesh_ptr_->getVertex(vm_s);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);

				shiftcoordv1(_x_n, _y_n, v_src_min);
				shiftcoordv4(_x_n, _y_n, v_src_min4);

				if (bedge_s != -1)
				{
					shiftbedgecoordv1(bedge_s, v_src_min);
					shiftbedgecoordv4(bedge_s, v_src_min4);
				}
			}
			//fill
			if ((centerv3(v_fld_plu3) - centerv3(v_src_ct3)).Norm()<1e-5)
			{
				Zpp = MDZppSingular(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
				Zpp = MDZppkernel(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);


			if ((centerv3(v_fld_min3) - centerv3(v_src_ct3)).Norm()<1e-5)
			{
				Zmp = MDZmpSingular(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
				Zmp = MDZmpkernel(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);

			if (vm_s != -1)
			{
				Zpm = MDZpmkernel(v_fld_plu3, v_src_min4, v_fld_plu, v_src_min, v_src_ct3, vm_s, epsp, epsm);
				Zmm = MDZmmkernel(v_fld_min3, v_src_min4, v_fld_min, v_src_min, v_src_ct3, vm_s, epsp, epsm);
			}

			Z_coup(f + end_tf - start_tf - start_ef, s - start_ts) = coef * l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f + end_tf - start_tf - start_ef, s - start_ts) = coef * l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void SED_CM_FMM::fillCOUPZMM_MP_B(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup)
{
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	Triangle tri_plu, tri_min;
	value_t l_fld, l_src;//场与源的公共边长度

	int sedid_f, sedid_s;
	int bedge_f, bedge_s;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;
	int arrID_f, arrID_s;
	arrID_f = _x_m * Array_y + _y_m;
	arrID_s = _x_n * Array_y + _y_n;

	ct_ptr_->getSEDCommonTriangleSize(_k_m, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(_k_m, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(_k_n, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(_k_n, start_es, end_es);


	Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm;
	Complex coef = (J0 * k_ * Z0) / PI4;

	//填充ZMM矩阵
	//tool::BarAndPercent bar_perc4;
	//coef = (J0*k_*Z0) / PI4;
	//#pragma omp parallel for
	for (int f = start_ef; f < end_ef; f++)
	{
		//bar_perc4(f - start_ef + 1, end_ef - start_ef);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld, sedid_f, bedge_f);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv1(_x_m, _y_m, v_fld_min);
		shiftcoordv3(_x_m, _y_m, v_fld_plu3);
		shiftcoordv3(_x_m, _y_m, v_fld_min3);

		if (bedge_f != -1)
		{
			shiftbedgecoordv1(bedge_f, v_fld_min);
			shiftbedgecoordv3(bedge_f, v_fld_min3);
		}
		for (int s = start_es + com_ce; s < end_es; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src, sedid_s, bedge_s);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv1(_x_n, _y_n, v_src_min);
			shiftcoordv3(_x_n, _y_n, v_src_plu3);
			shiftcoordv3(_x_n, _y_n, v_src_min3);

			if (bedge_s != -1)
			{
				shiftbedgecoordv1(bedge_s, v_src_min);
				shiftbedgecoordv3(bedge_s, v_src_min3);
			}
			//field face+ <-->  source face+
			if ((f_fld_plu == f_src_plu) && (arrID_f == arrID_s))
			{
				Zpp = MMZppSingular(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			else
			{
				Zpp = MMZppKernel(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			//field face+ <--> source face-
			if (f_fld_plu == f_src_min)
			{
				if ((centerv3(v_fld_plu3) - centerv3(v_src_min3)).Norm() < 1e-5)
				{
					Zpm = MMZpmSingular(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
				}
				else
				{
					Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
				}
			}
			else
			{
				Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
			}
			//field face- <--> source face+
			if (f_fld_min == f_src_plu)
			{
				if ((centerv3(v_fld_min3) - centerv3(v_src_plu3)).Norm() < 1e-5)
				{
					Zmp = MMZmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
				}
				else
				{
					Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
				}
			}
			else
			{
				Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
			}
			//field face- <--> source face-
			if (f_fld_min == f_src_min)
			{
				if ((centerv3(v_fld_min3) - centerv3(v_src_min3)).Norm() < 1e-5)
				{
					Zmm = MMZmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
				}
				else
				{
					Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
				}
			}
			else
			{
				Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
			}
			Z_coup(f + end_tf - start_tf - start_ef, s + end_ts - start_ts - start_es) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f + end_tf - start_tf - start_ef, s + end_ts - start_ts - start_es) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
	for (int f = start_ef + com_ce; f < end_ef; f++)
	{
		//bar_perc4(f - start_ef + 1, end_ef - start_ef);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld, sedid_f, bedge_f);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv1(_x_m, _y_m, v_fld_min);
		shiftcoordv3(_x_m, _y_m, v_fld_plu3);
		shiftcoordv3(_x_m, _y_m, v_fld_min3);

		if (bedge_f != -1)
		{
			shiftbedgecoordv1(bedge_f, v_fld_min);
			shiftbedgecoordv3(bedge_f, v_fld_min3);
		}
		for (int s = start_es; s < end_es; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src, sedid_s, bedge_s);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv1(_x_n, _y_n, v_src_min);
			shiftcoordv3(_x_n, _y_n, v_src_plu3);
			shiftcoordv3(_x_n, _y_n, v_src_min3);

			if (bedge_s != -1)
			{
				shiftbedgecoordv1(bedge_s, v_src_min);
				shiftbedgecoordv3(bedge_s, v_src_min3);
			}
			//field face+ <-->  source face+
			if ((f_fld_plu == f_src_plu) && (arrID_f == arrID_s))
			{
				Zpp = MMZppSingular(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			else
			{
				Zpp = MMZppKernel(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			//field face+ <--> source face-
			if (f_fld_plu == f_src_min)
			{
				if ((centerv3(v_fld_plu3) - centerv3(v_src_min3)).Norm() < 1e-5)
				{
					Zpm = MMZpmSingular(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
				}
				else
				{
					Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
				}
			}
			else
			{
				Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
			}
			//field face- <--> source face+
			if (f_fld_min == f_src_plu)
			{
				if ((centerv3(v_fld_min3) - centerv3(v_src_plu3)).Norm() < 1e-5)
				{
					Zmp = MMZmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
				}
				else
				{
					Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
				}
			}
			else
			{
				Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
			}
			//field face- <--> source face-
			if (f_fld_min == f_src_min)
			{
				if ((centerv3(v_fld_min3) - centerv3(v_src_min3)).Norm() < 1e-5)
				{
					Zmm = MMZmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
				}
				else
				{
					Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
				}
			}
			else
			{
				Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
			}
			Z_coup(f + end_tf - start_tf - start_ef, s + end_ts - start_ts - start_es) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f + end_tf - start_tf - start_ef, s + end_ts - start_ts - start_es) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
	for (int f = start_ef + com_ce; f < end_ef; f++)
	{
		//bar_perc4(f - start_ef + 1, end_ef - start_ef);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld, sedid_f, bedge_f);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);

		shiftcoordv1(_x_m, _y_m, v_fld_plu);
		shiftcoordv1(_x_m, _y_m, v_fld_min);
		shiftcoordv3(_x_m, _y_m, v_fld_plu3);
		shiftcoordv3(_x_m, _y_m, v_fld_min3);

		if (bedge_f != -1)
		{
			shiftbedgecoordv1(bedge_f, v_fld_min);
			shiftbedgecoordv3(bedge_f, v_fld_min3);
		}
		for (int s = start_es + com_ce; s < end_es; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src, sedid_s, bedge_s);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);

			shiftcoordv1(_x_n, _y_n, v_src_plu);
			shiftcoordv1(_x_n, _y_n, v_src_min);
			shiftcoordv3(_x_n, _y_n, v_src_plu3);
			shiftcoordv3(_x_n, _y_n, v_src_min3);

			if (bedge_s != -1)
			{
				shiftbedgecoordv1(bedge_s, v_src_min);
				shiftbedgecoordv3(bedge_s, v_src_min3);
			}
			//field face+ <-->  source face+
			if ((f_fld_plu == f_src_plu) && (arrID_f == arrID_s))
			{
				Zpp = MMZppSingular(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			else
			{
				Zpp = MMZppKernel(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			//field face+ <--> source face-
			if (f_fld_plu == f_src_min)
			{
				if ((centerv3(v_fld_plu3) - centerv3(v_src_min3)).Norm() < 1e-5)
				{
					Zpm = MMZpmSingular(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
				}
				else
				{
					Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
				}
			}
			else
			{
				Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
			}
			//field face- <--> source face+
			if (f_fld_min == f_src_plu)
			{
				if ((centerv3(v_fld_min3) - centerv3(v_src_plu3)).Norm() < 1e-5)
				{
					Zmp = MMZmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
				}
				else
				{
					Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
				}
			}
			else
			{
				Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
			}
			//field face- <--> source face-
			if (f_fld_min == f_src_min)
			{
				if ((centerv3(v_fld_min3) - centerv3(v_src_min3)).Norm() < 1e-5)
				{
					Zmm = MMZmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
				}
				else
				{
					Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
				}
			}
			else
			{
				Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
			}
			Z_coup(f + end_tf - start_tf - start_ef, s + end_ts - start_ts - start_es) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
			//Z_coup_eigen(f + end_tf - start_tf - start_ef, s + end_ts - start_ts - start_es) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void SED_CM_FMM::fillSEDV()
{
	int start_t, end_t, start_e, end_e;
	int nvp, nvm, fp, fm;
	int sedid, bedge;
	VectorR3 vp, vm, vp4[4], vm4[4];
	VectorR3 vp3[3], vm3[3];
	value_t area, ln;
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	int pol_from, pol_delta = 90, pol_to;
	if (multiInc.polarization == "Theta")
	{
		pol_from = 0;
		pol_to = 0;
	}
	else if (multiInc.polarization == "Phi")
	{
		pol_from = 90;
		pol_to = 90;
	}
	else
	{
		pol_from = 0;
		pol_to = 90;
	}
	tool::BarAndPercent bar_perc1;
	int theta_num = (multiInc.PW_theta_to - multiInc.PW_theta_from) / multiInc.PW_theta_delta + 1;
	int phi_num = (multiInc.PW_phi_to - multiInc.PW_phi_from) / multiInc.PW_phi_delta + 1;
	int k = 0;
	for (int pol = pol_from; pol <= pol_to; pol += pol_delta, k++)
	{
		int i = 0;
		for (value_t t = multiInc.PW_theta_from; t <= multiInc.PW_theta_to; t += multiInc.PW_theta_delta, i++)
		{
			bar_perc1(i + 1, theta_num);
			value_t theta = t * DegreesToRadians;
			int j = 0;
			for (value_t p = multiInc.PW_phi_from; p <= multiInc.PW_phi_to; p += multiInc.PW_phi_delta, j++)
			{
				value_t phi = p * DegreesToRadians;
				PW_inc_k = VectorR3(-sin(theta) * cos(phi), -sin(theta) * sin(phi), -cos(theta));
				PW_inc_e = VectorR3(-cos(pol)*cos(theta)*cos(phi) - sin(pol)*sin(phi), -cos(pol)*cos(theta)*sin(phi) + sin(pol)*cos(phi), cos(pol)*sin(theta));
				for (int id = 0; id < 9; id++)
				{
					ct_ptr_->getSEDCommonTriangleSize(id, start_t, end_t);
					ce_ptr_->getSEDCommonEdgeSize(id, start_e, end_e);
					for (int m = start_t; m < end_t; m++)
					{
						ct_ptr_->getCommonTriangle(m, nvp, nvm, fp, fm, area, sedid, bedge);
						vp = mesh_ptr_->getVertex(nvp);
						tet_plu = mesh_ptr_->getTetrahedronRef(fp);
						tet_plu.GetVertices(vp4[0], vp4[1], vp4[2], vp4[3]);

						shiftcoordv1(sedid, vp);
						shiftcoordv4(sedid, vp4);
						if (nvm != -1)
						{
							vm = mesh_ptr_->getVertex(nvm);
							tet_min = mesh_ptr_->getTetrahedronRef(fm);
							tet_min.GetVertices(vm4[0], vm4[1], vm4[2], vm4[3]);

							shiftcoordv1(sedid, vm);
							shiftcoordv4(sedid, vm4);

							if (bedge != -1)
							{
								shiftbedgecoordv1(bedge, vm);
								shiftbedgecoordv4(bedge, vm4);
							}
						}
						Complex vk = SEDDVKernel(vp4, vm4, vp, vm, nvm, PW_inc_k, PW_inc_e);
						V_sed(m + start_e, k*theta_num*phi_num + i * phi_num + j) = area * vk / 3.0f;
						V_sed_eigen(m + start_e, k*theta_num*phi_num + i * phi_num + j) = area * vk / 3.0f;
					}

					////////////////**************************************************///////////////////

					for (int n = start_e; n < end_e; n++)
					{
						ce_ptr_->getCommonEdge(n, nvp, nvm, fp, fm, ln, sedid, bedge);
						vp = mesh_ptr_->getVertex(nvp);
						tri_plu = mesh_ptr_->getTriangleRef(fp);
						tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);

						vm = mesh_ptr_->getVertex(nvm);
						tri_min = mesh_ptr_->getTriangleRef(fm);
						tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

						/*shiftcoordv1(sedid, vp);
						shiftcoordv1(sedid, vm);
						shiftcoordv3(sedid, vp3);
						shiftcoordv3(sedid, vm3);*/

						shiftcoordv1(sedid, vp, scale_x);
						shiftcoordv1(sedid, vm, scale_x);
						shiftcoordv3(sedid, vp3, scale_x);
						shiftcoordv3(sedid, vm3, scale_x);
						ln *= pow(scale_x, sedid / 3);

						if (bedge != -1)
						{
							shiftbedgecoordv1(bedge, vm);
							shiftbedgecoordv3(bedge, vm3);
						}

						Complex vk = SEDMVKernel(vp3, vm3, vp, vm, PW_inc_k, PW_inc_e);
						V_sed(n + end_t, k*theta_num*phi_num + i * phi_num + j) = 0.5f * ln * vk;
						V_sed_eigen(n + end_t, k*theta_num*phi_num + i * phi_num + j) = 0.5f * ln * vk;
					}
				}
			}


		}
	}

}

void SED_CM_FMM::fillSEDV_rad()
{
	int nv1, nv2, sedid;

	tool::BarAndPercent bar_perc;
	if (ex_num == 1)
	{
		for (int u = 0; u < unknowns_e; ++u)
		{
			bar_perc(u + 1, unknowns_e);

			ce_ptr_->getCommonEdge(u, nv1, nv2, sedid);
			auto length = ce_ptr_->getCommonEdgeLength(u);

			for (const auto& elem : exc_edge_)
			{
				if ((elem.first == nv1 && elem.second == nv2) || (elem.first == nv2 && elem.second == nv1))
				{
					//Qcout << "匹配到激励边！！！艹！" << std::endl;
					value_t phase;
					phase = (Phase_0 + (sedid / 3)*Phase_x + (sedid % 3)*Phase_y)*DegreesToRadians;
					if (elem.first == nv1)
						V_sed_eigen(u + unknowns_t, multiInc.PW_num) = length;// *exp(J0*phase);
					else
						V_sed_eigen(u + unknowns_t, multiInc.PW_num) = -length;// *exp(J0*phase);
																			   //break;
				}
			}
		}
	}
	else
	{
		for (int e = 0; e < exc_edge_.size(); e++)
		{
			int cur = 0;
			for (int u = 0; u < unknowns_e; ++u)
			{
				bar_perc(u + 1, unknowns_e);

				ce_ptr_->getCommonEdge(u, nv1, nv2, sedid);
				auto length = ce_ptr_->getCommonEdgeLength(u);

				if ((exc_edge_[e].first == nv1 && exc_edge_[e].second == nv2) || (exc_edge_[e].first == nv2 && exc_edge_[e].second == nv1))
				{
					//Qcout << "匹配到激励边！！！艹！" << std::endl;
					value_t phase;
					phase = (Phase_0 + (sedid / 3)*Phase_x + (sedid % 3)*Phase_y)*DegreesToRadians;

					if (exc_edge_[e].first == nv1)
						V_sed_eigen(u + unknowns_t, multiInc.PW_num + cur) = length;// *exp(J0*phase);
					else
						V_sed_eigen(u + unknowns_t, multiInc.PW_num + cur) = -length;// *exp(J0*phase);

					cur++;
					continue;
				}
			}
		}
	}


}

Complex SED_CM_FMM::SEDDVKernel(VectorR3 * vp4, VectorR3 * vm4, VectorR3 & vp, VectorR3 & vm, int & _nvm, VectorR3 &inc_k, VectorR3 &inc_e)
{
	VectorR3 vgp[4], vgm[4];
	VectorR3 rou_p, rou_m; //ρ+-
	Complex Ve;
	value_t kp = 0.0f, km = 0.0f;
	//kp = (epsp - 1.0f) / epsp;
	//km = (epsm - 1.0f) / epsm;

	Gauss4Point(vp4[0], vp4[1], vp4[2], vp4[3], vgp);
	Gauss4Point(vm4[0], vm4[1], vm4[2], vm4[3], vgm);

	VectorC3 Ei(0, 0, 0);
	Complex Vep(0, 0), Vem(0, 0), G(0, 0);

	for (int a = 0; a < 4; a++)
	{
		//////////////////////////////////////////////////////////////////////////
		rou_p = vgp[a] - vp;
		G = exp(-J0 * k_ * (vgp[a] ^ inc_k));

		Ei = G * inc_e;
		Vep += wt4[a] * (rou_p ^ Ei);

		//////////////////////////////////////////////////////////////////////////
		if (_nvm != -1)
		{
			rou_m = vm - vgm[a];
			G = exp(-J0 * k_ * (vgm[a] ^ inc_k));

			Ei = G * inc_e;
			Vem += wt4[a] * (rou_m ^ Ei);
		}

	}

	Ve = Vep + Vem;

	return Ve;
}

Complex SED_CM_FMM::SEDMVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &vp, VectorR3 &vm, VectorR3 &inc_k, VectorR3 &inc_e)
{
	VectorR3 vgp[7], vgm[7];
	VectorR3 rou_p, rou_m; //ρ+-
	Complex Ve;

	Gauss7Point(vp3[0], vp3[1], vp3[2], vgp);
	Gauss7Point(vm3[0], vm3[1], vm3[2], vgm);

	VectorC3 Ei(0, 0, 0)/*, Hi(0, 0, 0)*/;
	Complex Vep(0, 0), Vem(0, 0)/*, Vmp(0, 0), Vmm(0, 0)*/, G(0, 0);

	for (int a = 0; a < 7; a++)
	{
		//////////////////////////////////////////////////////////////////////////
		rou_p = vgp[a] - vp;
		G = exp(-J0 * k_ * (vgp[a] ^ inc_k));

		Ei = G * inc_e;
		Vep += w7[a] * (rou_p ^ Ei);

		//////////////////////////////////////////////////////////////////////////
		rou_m = vm - vgm[a];
		G = exp(-J0 * k_ * (vgm[a] ^ inc_k));

		Ei = G * inc_e;
		Vem += w7[a] * (rou_m ^ Ei);
	}

	Ve = (Vep + Vem);

	return Ve;
	//后面要乘以 （0.5*L）
}

void SED_CM_FMM::solveIsed()
{
	/*Qofstream outputZ(dir_ + "/matrix_Z.txt");
	Z_sed.save(outputZ, arma::arma_ascii);
	outputZ.flush();
	Qofstream outputV(dir_ + "/matrix_V.txt");
	V_sed.save(outputV, arma::arma_ascii);
	outputV.flush();*/
	//Qcx_vec I_temp, V_temp;
	CVector I_temp, V_temp;
	//I_temp.resize(unknowns_);
	Qcout << multiInc.PW_num << std::endl;

	int useless = 0;
	tool::BarAndPercent bar;
	/*for (int i = 0; i < multiInc.PW_num; i++)
	{
	bar(i + 1, multiInc.PW_num);
	//I_temp = I_sed.col(i);
	//if (!arma::solve(I_temp,Z_sed, V_sed.col(i),arma::solve_opts::fast))
	//throw std::runtime_error("fail solving matrix equation");
	//I_sed.col(i) = I_temp;

	V_temp = V_sed_eigen.col(i);
	I_temp = Z_sed_eigen.lu().solve(V_temp);
	I_sed_eigen.col(i) = I_temp;

	}*/
	//Eigen::BiCGSTAB<CMatrix> bicg;
	//bicg.compute(Z_sed_eigen);
	//I_sed_eigen = bicg.solve(V_sed_eigen);
	I_sed_eigen = Z_sed_eigen.partialPivLu().solve(V_sed_eigen);
}

void SED_CM_FMM::fillREDZ()
{
	int id_m, id_n;
	int start_t, end_t, start_e, end_e;
	int p = 0;
	Z_red.set_size(Array_x*Array_y*svd_k, Array_x*Array_y*svd_k);
	Z_red_eigen.resize(Array_x*Array_y*svd_k, Array_x*Array_y*svd_k);
	tool::BarAndPercent bar_perc;

	for (int x_m = 0; x_m < Array_x; x_m++)
	{
		for (int y_m = 0; y_m < Array_y; y_m++)
		{
			id_m = ArrayToSED(x_m, y_m);
			for (int x_n = 0; x_n < Array_x; x_n++)
			{
				for (int y_n = 0; y_n < Array_y; y_n++)
				{
					id_n = ArrayToSED(x_n, y_n);
					REDZKernel(x_m, y_m, x_n, y_n);// create coupling matrix "Z_coup" between two blocks

					for (int k_m = 0; k_m < svd_k; k_m++)//create the reduction matrix
					{
						auto& I_m = I_svd_eigen[id_m].col(k_m);
						int f = (x_m*Array_y + y_m)*svd_k + k_m;

						for (int k_n = 0; k_n < svd_k; k_n++)
						{
							auto& I_n = I_svd_eigen[id_n].col(k_n);
							int s = (x_n*Array_y + y_n)*svd_k + k_n;

							//Qcx_mat z_ele = I_m.t()*Z_coup*I_n;
							CMatrix z_ele = I_m.transpose()*Z_coup_eigen*I_n;
							//Z_red(f, s, arma::size(z_ele)) = z_ele;
							Z_red_eigen.block(f, s, z_ele.rows(), z_ele.cols()) = z_ele;
							bar_perc(p + 1, Array_x*Array_y*svd_k*Array_x*Array_y*svd_k);
							p = p + 1;
						}
					}

				}
			}


		}
	}
}

void SED_CM_FMM::fillREDV()
{
	V_red.set_size(Array_x*Array_y*svd_k);
	V_red_eigen.resize(Array_x*Array_y*svd_k, 1);
	for (int x = 0; x < Array_x; x++)
	{
		for (int y = 0; y < Array_y; y++)
		{
			int id = ArrayToSED(x, y);
			REDVKernel(x, y, id);

			for (int k = 0; k < svd_k; k++)
			{
				int f = (x*Array_y + y)*svd_k + k;
				auto& I_m = I_svd_eigen[id].col(k);
				auto& V_ele = I_m.transpose()*V_temp_eigen;
				V_red_eigen.block(f, 0, V_ele.rows(), V_ele.cols()) = V_ele;
			}
		}
	}
}

void SED_CM_FMM::fillREDV_MP()
{
	auto task_num = Task_assign2.size() - 1;
	V_red_eigen.resize(Array_x*Array_y*svd_k, 1);

	for (int t = 0; t < task_num; t++)
	{
		int begin = Task_assign3[t];
		int end = Task_assign3[t + 1];

		if (begin < end)
		{
			auto task = std::bind(&SED_CM_FMM::REDV_MP_Kernel, this, begin, end);
			pool_.submit(task);
		}
	}
	Qcout << "Finish V task submit" << std::endl;
	pool_.interactiveRun();

}

void SED_CM_FMM::REDV_MP_Kernel(int begin, int end)
{
	for (int u = begin; u < end; u++)
	{
		int x = u / Array_y;
		int y = u % Array_y;

		int id = ArrayToSED(x, y);
		//REDVKernel_rad(x, y, id);
		REDVKernel(x, y, id);

		for (int k = 0; k < svd_k; k++)
		{
			int f = (x*Array_y + y)*svd_k + k;
			auto& I_m = I_svd_eigen[id].col(k);
			auto& V_ele = I_m.transpose()*V_temp_eigen;
			V_red_eigen.block(f, 0, V_ele.rows(), V_ele.cols()) = V_ele;
		}
	}
}

void SED_CM_FMM::fillREDV_rad()
{
	V_red.set_size(Array_x*Array_y*svd_k);
	V_red_eigen.resize(Array_x*Array_y*svd_k, 1);
	for (int x = 0; x < Array_x; x++)
	{
		for (int y = 0; y < Array_y; y++)
		{
			int id = ArrayToSED(x, y);
			REDVKernel_rad(x, y, id);

			for (int k = 0; k < svd_k; k++)
			{
				int f = (x*Array_y + y)*svd_k + k;
				auto& I_m = I_svd_eigen[id].col(k);
				auto& V_ele = I_m.transpose()*V_temp_eigen;
				V_red_eigen.block(f, 0, V_ele.rows(), V_ele.cols()) = V_ele;
			}
		}
	}
}

void SED_CM_FMM::REDZKernel(int &_x_m, int &_y_m, int &_x_n, int &_y_n)
{
	int k_f, k_s;
	int bedge_t, sedid_t, bedge_e, sedid_e;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;

	k_f = ArrayToSED(_x_m, _y_m);
	k_s = ArrayToSED(_x_n, _y_n);
	ct_ptr_->getSEDCommonTriangleSize(k_f, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(k_f, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(k_s, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(k_s, start_es, end_es);
	Z_coup.resize(end_tf - start_tf + end_ef - start_ef, end_ts - start_ts + end_es - start_es);
	Z_coup_eigen.resize(end_tf - start_tf + end_ef - start_ef, end_ts - start_ts + end_es - start_es);


	fillCOUPZDD(_x_m, _y_m, _x_n, _y_n, k_f, k_s);
	fillCOUPZDM(_x_m, _y_m, _x_n, _y_n, k_f, k_s);
	fillCOUPZMD(_x_m, _y_m, _x_n, _y_n, k_f, k_s);
	fillCOUPZMM(_x_m, _y_m, _x_n, _y_n, k_f, k_s);
}

CMatrix SED_CM_FMM::coupZ_MP(int &_x_m, int &_y_m, int &_x_n, int &_y_n)
{
	int k_f, k_s;
	int bedge_t, sedid_t, bedge_e, sedid_e;
	int start_tf, end_tf, start_ts, end_ts;
	int start_ef, end_ef, start_es, end_es;

	k_f = 0;
	k_s = 0;
	ct_ptr_->getSEDCommonTriangleSize(k_f, start_tf, end_tf);
	ce_ptr_->getSEDCommonEdgeSize(k_f, start_ef, end_ef);
	ct_ptr_->getSEDCommonTriangleSize(k_s, start_ts, end_ts);
	ce_ptr_->getSEDCommonEdgeSize(k_s, start_es, end_es);
	CMatrix Z_coup;
	Z_coup.resize(end_tf - start_tf + end_ef - start_ef, end_ts - start_ts + end_es - start_es);
	//Z_coup_eigen.resize(end_tf - start_tf + end_ef - start_ef, end_ts - start_ts + end_es - start_es);

	fillCOUPZDD_MP(_x_m, _y_m, _x_n, _y_n, k_f, k_s, Z_coup);
	fillCOUPZDM_MP(_x_m, _y_m, _x_n, _y_n, k_f, k_s, Z_coup);
	fillCOUPZMD_MP(_x_m, _y_m, _x_n, _y_n, k_f, k_s, Z_coup);
	fillCOUPZMM_MP(_x_m, _y_m, _x_n, _y_n, k_f, k_s, Z_coup);

	return Z_coup;
}

CMatrix SED_CM_FMM::coupZ_MP_COM(int &_x_m, int &_y_m, int &_x_n, int &_y_n)
{
	int k_f, k_s;

	k_f = 0;
	k_s = 0;

	CMatrix Z_coup_com;
	Z_coup_com.resize(com_ct + com_ce, com_ct + com_ce);
	//Z_coup_eigen.resize(end_tf - start_tf + end_ef - start_ef, end_ts - start_ts + end_es - start_es);

	fillCOUPZDD_MP_C(_x_m, _y_m, _x_n, _y_n, k_f, k_s, Z_coup_com);
	fillCOUPZDM_MP_C(_x_m, _y_m, _x_n, _y_n, k_f, k_s, Z_coup_com);
	fillCOUPZMD_MP_C(_x_m, _y_m, _x_n, _y_n, k_f, k_s, Z_coup_com);
	fillCOUPZMM_MP_C(_x_m, _y_m, _x_n, _y_n, k_f, k_s, Z_coup_com);

	return Z_coup_com;
}

void SED_CM_FMM::REDVKernel(int &_x, int &_y, int &_id)
{
	int nvp, nvm, fp, fm, bedge, sedid;
	VectorR3 vp, vm, vp4[4], vm4[4];
	VectorR3 vp3[3], vm3[3];
	value_t area, ln;
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;

	int start_t, end_t, start_e, end_e;

	ct_ptr_->getSEDCommonTriangleSize(_id, start_t, end_t);
	ce_ptr_->getSEDCommonEdgeSize(_id, start_e, end_e);

	V_temp.resize(end_t - start_t + end_e - start_e);
	V_temp_eigen.resize(end_t - start_t + end_e - start_e);

	for (int m = start_t; m < end_t; m++)
	{
		ct_ptr_->getCommonTriangle(m, nvp, nvm, fp, fm, area, sedid, bedge);

		vp = mesh_ptr_->getVertex(nvp);
		tet_plu = mesh_ptr_->getTetrahedronRef(fp);
		tet_plu.GetVertices(vp4[0], vp4[1], vp4[2], vp4[3]);

		shiftcoordv1(_x, _y, vp);
		shiftcoordv4(_x, _y, vp4);
		if (nvm != -1)
		{
			vm = mesh_ptr_->getVertex(nvm);
			tet_min = mesh_ptr_->getTetrahedronRef(fm);
			tet_min.GetVertices(vm4[0], vm4[1], vm4[2], vm4[3]);

			shiftcoordv1(_x, _y, vm);
			shiftcoordv4(_x, _y, vm4);

			if (bedge != -1)
			{
				shiftbedgecoordv1(bedge, vm);
				shiftbedgecoordv4(bedge, vm4);
			}
		}

		Complex vk = DVKernel(vp4, vm4, vp, vm, nvm);
		V_temp(m - start_t) = area * vk / 3.0f;
		V_temp_eigen(m - start_t) = area * vk / 3.0f;
	}

	////////////////**************************************************///////////////////

	for (int m = start_e; m < end_e; m++)
	{
		ce_ptr_->getCommonEdge(m, nvp, nvm, fp, fm, ln, sedid, bedge);
		vp = mesh_ptr_->getVertex(nvp);
		tri_plu = mesh_ptr_->getTriangleRef(fp);
		tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);

		vm = mesh_ptr_->getVertex(nvm);
		tri_min = mesh_ptr_->getTriangleRef(fm);
		tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

		shiftcoordv1(_x, _y, vp, scale_x);
		shiftcoordv1(_x, _y, vm, scale_x);
		shiftcoordv3(_x, _y, vp3, scale_x);
		shiftcoordv3(_x, _y, vm3, scale_x);
		ln *= pow(scale_x, _x);

		if (bedge != -1)
		{
			shiftbedgecoordv1(bedge, vm);
			shiftbedgecoordv3(bedge, vm3);
		}

		Complex vk = MVKernel(vp3, vm3, vp, vm);
		V_temp(m - start_e + end_t - start_t) = 0.5f * ln * vk;
		V_temp_eigen(m - start_e + end_t - start_t) = 0.5f * ln * vk;
	}
}

void SED_CM_FMM::REDVKernel_rad(int &_x, int &_y, int &_id)
{
	int nvp, nvm, fp, fm, bedge, sedid, nv1, nv2;
	VectorR3 vp, vm, vp4[4], vm4[4];
	VectorR3 vp3[3], vm3[3];
	value_t area, ln;
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;

	int start_t, end_t, start_e, end_e;

	ct_ptr_->getSEDCommonTriangleSize(_id, start_t, end_t);
	ce_ptr_->getSEDCommonEdgeSize(_id, start_e, end_e);

	V_temp.resize(end_t - start_t + end_e - start_e);
	V_temp_eigen.setZero(end_t - start_t + end_e - start_e);

	/*for (int m = start_t; m < end_t; m++)
	{
	ct_ptr_->getCommonTriangle(m, nvp, nvm, fp, fm, area, sedid, bedge);

	vp = mesh_ptr_->getVertex(nvp);
	tet_plu = mesh_ptr_->getTetrahedronRef(fp);
	tet_plu.GetVertices(vp4[0], vp4[1], vp4[2], vp4[3]);

	shiftcoordv1(_x, _y, vp);
	shiftcoordv4(_x, _y, vp4);
	if (nvm != -1)
	{
	vm = mesh_ptr_->getVertex(nvm);
	tet_min = mesh_ptr_->getTetrahedronRef(fm);
	tet_min.GetVertices(vm4[0], vm4[1], vm4[2], vm4[3]);

	shiftcoordv1(_x, _y, vm);
	shiftcoordv4(_x, _y, vm4);

	if (bedge != -1)
	{
	shiftbedgecoordv1(bedge, vm);
	shiftbedgecoordv4(bedge, vm4);
	}
	}

	Complex vk = DVKernel(vp4, vm4, vp, vm, nvm);
	V_temp(m - start_t) = area * vk / 3.0f;
	V_temp_eigen(m - start_t) = area * vk / 3.0f;
	}*/

	////////////////**************************************************///////////////////

	for (int m = start_e; m < end_e; m++)
	{
		ce_ptr_->getCommonEdge(m, nv1, nv2);

		auto length = ce_ptr_->getCommonEdgeLength(m);

		for (const auto& elem : exc_edge_)
		{
			if ((elem.first == nv1 && elem.second == nv2) || (elem.first == nv2 && elem.second == nv1))
			{
				value_t mag;
				mag = 1.0f;
				value_t t = 1.0f;

				if (_x < 9.5)
					mag = _x * t + 1;
				else
					mag = (1 + 9 * t) - t * (_x - 10);


				value_t phase;
				if (_x == 2)
					phase = (Phase_0 + _x * Phase_x + _y * Phase_y)*DegreesToRadians;
				else
					phase = (Phase_0 + _x * Phase_x + _y * Phase_y)*DegreesToRadians;

				Complex phase_s;
				phase_s = exp(J0*phase);

				if (elem.first == nv1)
					V_temp_eigen(m - start_e + end_t - start_t) = length * phase_s*mag;
				else
					V_temp_eigen(m - start_e + end_t - start_t) = -length * phase_s*mag;
			}
		}
		//V_temp(m - start_e + end_t - start_t) = 0.5f * ln * vk;
		//V_temp_eigen(m - start_e + end_t - start_t) = 0.5f * ln * vk;
	}
}

void SED_CM_FMM::fillSVDI()
{
	//Qcx_mat U, V, I_pw;
	CMatrix U, V, I_pw;
	Qvec s;
	int sedid, N_e, N_pw, col_tot;
	int start_t, end_t, start_e, end_e;
	col_tot = I_sed_eigen.cols();
	int k = 0;
	value_t eps = 0;
	for (int id = 0; id < 9; id++)
	{
		ct_ptr_->getSEDCommonTriangleSize(id, start_t, end_t);
		ce_ptr_->getSEDCommonEdgeSize(id, start_e, end_e);
		N_e = end_t - start_t + end_e - start_e;

		if (multiInc.PW_num != 0)
		{
			if (rtype_ == policy::ResultType::Rad)//c
			{
				I_pw.resize(N_e, col_tot - ex_num);
				I_pw = I_sed_eigen.block(start_e, 0, N_e, col_tot - ex_num);
			}
			else
			{
				I_pw.resize(N_e, col_tot);
				I_pw = I_sed_eigen.block(start_e, 0, N_e, col_tot);
			}

			Eigen::JacobiSVD<CMatrix> svd(I_pw, Eigen::ComputeThinU);

			auto& svd_val = svd.singularValues();
			auto& svd_mat = svd.matrixU();
			auto max_val = svd_val.maxCoeff();

			//auto svd_matn = svd_mat * (svd_val.asDiagonal());//test

			while (1)
			{
				if ((svd_val[k] / max_val) < svd_threshold)
					break;

				k = k + 1;

				if (k == svd_val.size())
					break;
			}
			k = k + ex_num;//c
			svd_k = k >= svd_k ? k : svd_k;

			if (rtype_ == policy::ResultType::Rad)
			{
				CMatrix svd_mat_temp;
				svd_mat_temp.resize(svd_mat.rows(), svd_mat.cols() + ex_num);
				svd_mat_temp.block(0, 0, N_e, ex_num) = I_sed_eigen.block(start_e, col_tot - ex_num, N_e, ex_num);
				svd_mat_temp.rightCols(svd_mat.cols()) = svd_mat;
				I_svd_eigen.push_back(svd_mat_temp);
			}
			else
				I_svd_eigen.push_back(svd_mat);
				//I_svd_eigen.push_back(I_pw);//test
		}
		else
		{
			I_svd_eigen.push_back(I_sed_eigen.block(start_e, 0, N_e, col_tot));
			svd_k = col_tot;
		}

		Qcout << "the k post SVD is :" << '\n' << svd_k << std::endl;
		k = 0;//return to 0
	}

	//svd_k = 3;//test

}

void SED_CM_FMM::constructI()
{
	int id;
	int start_t, end_t, start_e, end_e;
	int blk_unk;
	int blk_num;
	//Qcx_vec ex_coff;
	for (int x = 0; x < Array_x; x++)
	{
		for (int y = 0; y < Array_y; y++)
		{
			id = ArrayToSED(x, y);
			ct_ptr_->getSEDCommonTriangleSize(id, start_t, end_t);
			ce_ptr_->getSEDCommonEdgeSize(id, start_e, end_e);

			blk_unk = end_t - start_t + end_e - start_e;
			blk_num = x * Array_y + y;

			//Qcx_vec ex_coff = I_red.subvec(blk_num*svd_k, (blk_num + 1)*svd_k-1);
			//Qcx_vec Iblk = (I_svd[id].cols(0, svd_k - 1))*ex_coff;
			//I_tot.push_back(Iblk);

			CVector Iblk_eigen = (I_svd_eigen[id].leftCols(svd_k))*I_red_eigen.segment(blk_num*svd_k, svd_k);
			//Iblk_eigen.segment(0, end_t - start_t) = (J0*omiga / Z0)*Iblk_eigen.segment(0, end_t - start_t);
			I_tot_eigen.push_back(Iblk_eigen);
		}
	}

}

void SED_CM_FMM::getJCM()
{
	CMatrix R = (Z_sed_eigen + Z_sed_eigen.conjugate()) / 2.0f;
	CMatrix X = (Z_sed_eigen - Z_sed_eigen.conjugate()) / (2.0f*J0);

	Eigen::FullPivLU<CMatrix> lu(X);
	CMatrix C = X.inverse()*R;

	Eigen::ComplexEigenSolver<CMatrix> Jcm;
	Jcm.compute(C);
	CMatrix Jn = Jcm.eigenvectors();
	CVector CMV = Jcm.eigenvalues().cwiseInverse();//(1/lamda)->lamda
	auto cm_tot = CMV.size();

	int k = 0;
	while (1)
	{
		//value_t ms0 = 1.0f/ (abs(abs(CMV[cm_tot-1].real())*J0+1.0f));
		//value_t msk = 1.0f / (abs(abs(CMV[cm_tot-1-k].real())*J0+1.0f));
		value_t msk = 1.0f / (abs(CMV[cm_tot - 1 - k] * J0 + 1.0f));
		if ((abs(msk)) < svd_threshold)
			break;
		k = k + 1;

		if (k == CMV.size())
			break;
	}
	svd_k = k >= svd_k ? k : svd_k;
	Qcout << "The number of Jcm is: " << svd_k << std::endl;
	int start_t, end_t, start_e, end_e;
	int N_e;
	CMatrix jcm;
	int col_jn = Jn.cols();

	for (int id = 0; id < 9; id++)
	{
		ct_ptr_->getSEDCommonTriangleSize(id, start_t, end_t);
		ce_ptr_->getSEDCommonEdgeSize(id, start_e, end_e);
		N_e = end_t - start_t + end_e - start_e;
		jcm.resize(N_e, svd_k);
		jcm = Jn.block(start_e, col_jn - svd_k - 1, N_e, svd_k);
		I_svd_eigen.push_back(jcm);
	}

	Qofstream outputJn(dir_ + "/martix_Jn.txt");
	outputJn << CMV;
	outputJn.flush();
	outputJn.close();

	/*Eigen::JacobiSVD<CMatrix> svd(R, Eigen::ComputeThinU | Eigen::ComputeThinV);
	auto& svd_val = svd.singularValues();
	auto& U = svd.matrixU();
	auto& V = svd.matrixV();
	auto max_val = svd_val.maxCoeff();

	int k=0;

	while (1)
	{
	if ((svd_val[k] / max_val) < svd_threshold)
	break;

	k = k + 1;

	if (k == svd_val.size())
	break;
	}

	k = k + ex_num;//c
	svd_k = k >= svd_k ? k : svd_k;
	//svd_k = 5;
	Qcout << "The number of Jcm is" << svd_k << std::endl;
	CMatrix A = U.transpose()*X*V;
	CMatrix A11 = A.topLeftCorner(svd_k, svd_k);
	CMatrix A12 = A.topRightCorner(svd_k, A.cols() - svd_k);
	CMatrix A22 = A.bottomRightCorner(A.rows() - svd_k, A.cols() - svd_k);
	//CMatrix W(svd_val.segment(0, svd_k).asDiagonal());
	CVector W11 = svd_val.segment(0, svd_k);
	CMatrix W11t(W11.cwiseInverse().cwiseSqrt().asDiagonal());//W11-1/2
	CMatrix A22i = A22.inverse();
	CMatrix C = W11t * (A11 - A12 * A22i*A12.transpose())*W11t;

	Eigen::ComplexEigenSolver<CMatrix> Jcm;
	Jcm.compute(C);
	auto & Y = Jcm.eigenvectors();

	CMatrix B1 = W11t*Y;
	CMatrix B2 = -A22i * A12.transpose()*B1;
	CMatrix B;
	B.resize(B1.rows() + B2.rows(), B1.cols());
	B.block(0, 0, B1.rows(), B1.cols()) = B1;
	B.block(B1.rows(), 0, B2.rows(), B2.cols()) = B2;
	CMatrix Jn = V * B;

	Qofstream outputJn(dir_ + "/martix_Jn.txt");
	outputJn << Z_sed_eigen;
	outputJn.flush();
	outputJn.close();

	int start_t, end_t, start_e, end_e;
	int N_e;
	CMatrix jcm;
	int col_jn = Jn.cols();

	for (int id = 0; id < 9; id++)
	{
	ct_ptr_->getSEDCommonTriangleSize(id, start_t, end_t);
	ce_ptr_->getSEDCommonEdgeSize(id, start_e, end_e);
	N_e = end_t - start_t + end_e - start_e;
	jcm.resize(N_e, col_jn);
	jcm = Jn.block(start_e, 0, N_e, col_jn);
	I_svd_eigen.push_back(jcm);
	}*/

	Qcout << "The CM current has been completed!" << std::endl;
	system("pause");

}

void SED_CM_FMM::buildblk()
{
	for (int idx = 0; idx < Array_x; idx++)
	{
		for (int idy = 0; idy < Array_y; idy++)
		{
			Blk blk(idx, idy);
			blk._id_sed = ArrayToSED(idx, idy);
			blk._id_near = ArrayToNear(idx, idy);
			BlockArray.emplace_back(blk);
		}
	}

	for (auto& pblk : BlockArray)
		findnearblk(pblk);

	for (int b = 1; b < BlockArray.size(); b++)
	{
		BlockArray[b].offset = svd_k * svd_k*(BlockArray[b - 1].near.size() + 1) + BlockArray[b - 1].offset;
		BlockArray[b].preoffset = svd_k * svd_k*(BlockArray[b - 1].near.size()) + BlockArray[b - 1].preoffset;
	}
		
}

void SED_CM_FMM::findnearblk(Blk& blk)
{
	int idx = blk._idx;
	int idy = blk._idy;
	int idm_near = blk._id_near;
	
	if (idx - 1 >= 0)
	{
		//(x-1,y)
		blk.near.emplace_back(&BlockArray[(idx - 1)*Array_y + idy]);

		if (idy + 1 < Array_y)//(x-1,y+1)
			blk.near.emplace_back(&BlockArray[(idx - 1)*Array_y + idy + 1]);
			
	}

	if (idy + 1 < Array_y)
	{
		//(x,y+1)
		blk.near.emplace_back(&BlockArray[idx*Array_y + idy + 1]);

		if (idx + 1 < Array_x)//(x+1,y+1)
			blk.near.emplace_back(&BlockArray[(idx + 1)*Array_y + idy + 1]);

	}

	if (idx + 1 < Array_x)
	{
		//(x+1,y)
		blk.near.emplace_back(&BlockArray[(idx + 1)*Array_y + idy]);

		if (idy - 1 >= 0)//(x+1,y-1)
			blk.near.emplace_back(&BlockArray[(idx + 1)*Array_y + idy - 1]);
	}

	if (idy - 1 >= 0)
	{
		//(x,y-1)
		blk.near.emplace_back(&BlockArray[idx*Array_y + idy - 1]);

		if (idx - 1 >= 0)//(x-1,y-1)
			blk.near.emplace_back(&BlockArray[(idx - 1)*Array_y + idy - 1]);
	}
}

void SED_CM_FMM::preCalculateKArray()
{
	gl_theta_.set(sam_theta_);
	gl_phi_.set(sam_phi_);

	std::vector<value_t> sint(sam_theta_), cost(sam_theta_);
	for (int i = 0; i < sam_theta_; ++i)
	{
		value_t theta = PIhalves * (gl_theta_.root(i) + 1);
		sint[i] = sin(theta);
		cost[i] = cos(theta);
	}

	std::vector<value_t> sinp(sam_phi_), cosp(sam_phi_);
	for (int i = 0; i < sam_phi_; ++i)
	{
		value_t phi = PI * (gl_phi_.root(i) + 1);
		sinp[i] = sin(phi);
		cosp[i] = cos(phi);
	}

	arr_k_.set_size(sam_theta_, sam_phi_);
	for (int p = 0; p < sam_phi_; ++p)
		for (int t = 0; t < sam_theta_; ++t)
			arr_k_(t, p) = VectorR3(sint[t] * cosp[p], sint[t] * sinp[p], cost[t]);
}

void SED_CM_FMM::preCalculateTrans()
{
	const size_t knum = sam_theta_ * sam_phi_;
	const Complex j_l[4] = { 1, -J0, -1, J0 };
	SphHankel sph_hankel;
	Legendre legendre;
	transfer_.resize(knum);
	for (size_t k = 0; k < knum; ++k)
	{
		transfer_[k].resize(2 * Array_x - 1, 2 * Array_y - 1);
		for (int xid = 1 - Array_x; xid < Array_x; ++xid)
		{
			for (int yid = 1 - Array_y; yid < Array_y; ++yid)
			{
				if (xid != 0 || yid != 0)
				{
					VectorR3 rvec(xid*Dx, yid*Dy, 0.0);
					const value_t vec_norm = rvec.Norm();
					const value_t kd_norm = k_ * vec_norm;
					const value_t kd_unit = (arr_k_(k) ^ rvec) / vec_norm;
					int minL = static_cast<int>(sam_theta_ < kd_norm ? sam_theta_ : kd_norm);
					Complex val(0, 0);
					for (int cl = 0; cl <= minL; ++cl)
						val += (2 * cl + 1.0f) * j_l[cl % 4] * sph_hankel(kd_norm, cl) * legendre(kd_unit, cl);
					transfer_[k](xid + Array_x - 1, yid + Array_y - 1) = val;
				}
			}
		}
	}
}

void SED_CM_FMM::preCalculateTrans_v2()
{
	const size_t knum = sam_theta_ * sam_phi_;
	transfer_.resize(knum);
	for(size_t k = 0; k < knum; ++k)
		transfer_[k].resize(2 * Array_x - 1, 2 * Array_y - 1);

	auto task_num = Task_assign2.size() - 1;
	for (int t = 0; t < task_num; t++)
	{
		int begin = Task_assign2[t];
		int end = Task_assign2[t + 1];
		if (begin < end)
		{
			auto task = std::bind(&SED_CM_FMM::CalculateTransKernel, this, begin, end);
			pool_.submit(task);
		}
	}
	pool_.interactiveRun();
}

void SED_CM_FMM::CalculateTransKernel(int begin, int end)
{
	const size_t knum = sam_theta_ * sam_phi_;
	const Complex j_l[4] = { 1, -J0, -1, J0 };
	SphHankel sph_hankel;
	Legendre legendre;
	
	for (int t = begin; t < end; t++)
	{
		int xid = SED_task2[t].first;
		int yid = SED_task2[t].second;

		for (size_t k = 0; k < knum; ++k)
		{
			if (xid != 0 || yid != 0)
			{
				VectorR3 rvec(xid*Dx, yid*Dy, 0.0);
				const value_t vec_norm = rvec.Norm();
				const value_t kd_norm = k_ * vec_norm;
				const value_t kd_unit = (arr_k_(k) ^ rvec) / vec_norm;
				int minL = static_cast<int>(sam_theta_ < kd_norm ? sam_theta_ : kd_norm);
				Complex val(0, 0);
				for (int cl = 0; cl <= minL; ++cl)
					val += (2 * cl + 1.0f) * j_l[cl % 4] * sph_hankel(kd_norm, cl) * legendre(kd_unit, cl);
				transfer_[k](xid + Array_x - 1, yid + Array_y - 1) = val;
			}
		}
	}
	
}

void SED_CM_FMM::preCalculateTrans(int& _x, int&_y,CVector& _transfer) const
{
	const size_t knum = sam_theta_ * sam_phi_;
	const Complex j_l[4] = { 1, -J0, -1, J0 };
	SphHankel sph_hankel;
	Legendre legendre;
	_transfer.resize(knum);
	//Qcout << "Begin calculate transfer..." << std::endl;
	VectorR3 rvec(Dx*_x, Dy*_y, 0);
	const value_t vec_norm = rvec.Norm();
	const value_t kd_norm = k_ * vec_norm;
	for (int k = 0; k < knum; k++)
	{
		const value_t kd_unit = (arr_k_(k) ^ rvec) / vec_norm;
		int minL = static_cast<int>(sam_theta_ < kd_norm ? sam_theta_ : kd_norm);
		Complex val(0, 0);
		for (int cl = 0; cl <= minL; ++cl)
			val += (2 * cl + 1.0f) * j_l[cl % 4] * sph_hankel(kd_norm, cl) * legendre(kd_unit, cl);
		_transfer(k) = val;
	}
	//Qcout << "End calculate transfer!!!" <<_transfer.size()<< std::endl;
}

void SED_CM_FMM::preCalculateRadAndRecv()
{
	/*const size_t box_num = Array_x;
	const size_t knum = sam_phi_ * sam_theta_;
	Qcout << knum << ' ' << com_ce << std::endl;
	int y = 0;
	tool::BarAndPercent bar_perc;   //
	for (int idx = 0; idx < 1; ++idx)
	{
		bar_perc(idx + 1, box_num);     //

		sed_rad_.set_size(com_ce, knum);
		sed_rec_.set_size(knum, com_ce);

		VectorR3 v_plu, v_min, nml_plu, nml_min;
		VectorR3 v_plu3[3], v_min3[3], vpg3[3], vmg3[3];

		for (size_t u = 0; u < com_ce; ++u)
		{
			bar_perc(u + 1, com_ce);
			//int u = i;
			int nvp, nvm, nfp, nfm;
			value_t length;
			ce_ptr_->getCommonEdge(u, nvp, nvm, nfp, nfm, length);
			v_plu = mesh_ptr_->getVertex(nvp);
			v_min = mesh_ptr_->getVertex(nvm);

			auto& tri_plu = mesh_ptr_->getTriangleRef(nfp);
			auto& tri_min = mesh_ptr_->getTriangleRef(nfm);
			tri_plu.getVertex(v_plu3[0], v_plu3[1], v_plu3[2]);
			tri_min.getVertex(v_min3[0], v_min3[1], v_min3[2]);
			nml_plu = tri_plu.getNormal();
			nml_min = tri_min.getNormal();

			Gauss3Point(v_plu3[0], v_plu3[1], v_plu3[2], vpg3);
			Gauss3Point(v_min3[0], v_min3[1], v_min3[2], vmg3);

			for (size_t k = 0; k < knum; ++k)
			{
				const auto tmp_rad = radiationFunction(vpg3, vmg3, v_plu, v_min, cen_box, arr_k_(k));
				//const auto tmp_rec = halfRecvFunction(vpg3, vmg3, v_plu, v_min, nml_plu, nml_min, center, arr_k_(k));
				//rad_mat(u, k) = (length / 2) * tmp_rad;
				//rec_mat(k, u) = (length / 2) * (conj(tmp_rad));
				sed_rad_(u, k) = (length / 2) * tmp_rad;
				sed_rec_(k, u) = (length / 2) * (conj(tmp_rad));
			}
		}
	}*/ // end for box
	const size_t knum = sam_phi_ * sam_theta_;
	for (int id = 0; id < 9; id++)
	{

		Qc3_mat rad_temp_, rec_temp_;
		rad_temp_.set_size(knum, svd_k);
		rec_temp_.set_size(svd_k, knum);

		for (int s_k = 0; s_k < svd_k; s_k++)
		{
			CVector I_temp = I_svd_eigen[id].col(s_k);

			for (int k = 0; k < knum; k++)
			{
				VectorC3 rad, rec;
				for (int u = 0; u < com_ce; u++)
				{
					int nvp, nvm, nfp, nfm;
					VectorR3 v_plu, v_min, nml_plu, nml_min;
					VectorR3 v_plu3[3], v_min3[3], vpg3[3], vmg3[3];
					value_t length;
					ce_ptr_->getCommonEdge(u, nvp, nvm, nfp, nfm, length);
					v_plu = mesh_ptr_->getVertex(nvp);
					v_min = mesh_ptr_->getVertex(nvm);

					auto& tri_plu = mesh_ptr_->getTriangleRef(nfp);
					auto& tri_min = mesh_ptr_->getTriangleRef(nfm);
					tri_plu.getVertex(v_plu3[0], v_plu3[1], v_plu3[2]);
					tri_min.getVertex(v_min3[0], v_min3[1], v_min3[2]);
					nml_plu = tri_plu.getNormal();
					nml_min = tri_min.getNormal();

					Gauss3Point(v_plu3[0], v_plu3[1], v_plu3[2], vpg3);
					Gauss3Point(v_min3[0], v_min3[1], v_min3[2], vmg3);

					const auto tmp_rad = radiationFunction(vpg3, vmg3, v_plu, v_min, cen_box, arr_k_(k));
					//const auto tmp_rec = halfRecvFunction(vpg3, vmg3, v_plu, v_min, nml_plu, nml_min, center, arr_k_(k));
					//rad_mat(u, k) = (length / 2) * tmp_rad;
					//rec_mat(k, u) = (length / 2) * (conj(tmp_rad));
					rad += (length / 2) * tmp_rad * I_temp(u);
					rec += (length / 2) * (conj(tmp_rad)) * I_temp(u);
				}
				rad_temp_(k, s_k) = rad;
				rec_temp_(s_k, k) = rec;      
			}
		}
		sed_rad_.push_back(rad_temp_);
		sed_rec_.push_back(rec_temp_);
	}

}

void SED_CM_FMM::preCalculateRadAndRecv_v2()
{
	const size_t knum = sam_phi_ * sam_theta_;
	Qc3_mat rad_com_, rec_com_;
	rad_com_.set_size(knum, com_ce);
	rec_com_.set_size(com_ce, knum);
	std::vector<VectorR3> dp;
	/*for (int k = 0; k < knum; k++)
	{
		for (int u = 0; u < com_ce; u++)
		{
			int nvp, nvm, nfp, nfm;
			VectorR3 v_plu, v_min, nml_plu, nml_min;
			VectorR3 v_plu3[3], v_min3[3], vpg3[3], vmg3[3];
			value_t length;
			ce_ptr_->getCommonEdge(u, nvp, nvm, nfp, nfm, length);
			v_plu = mesh_ptr_->getVertex(nvp);
			v_min = mesh_ptr_->getVertex(nvm);

			auto& tri_plu = mesh_ptr_->getTriangleRef(nfp);
			auto& tri_min = mesh_ptr_->getTriangleRef(nfm);
			tri_plu.getVertex(v_plu3[0], v_plu3[1], v_plu3[2]);
			tri_min.getVertex(v_min3[0], v_min3[1], v_min3[2]);
			nml_plu = tri_plu.getNormal();
			nml_min = tri_min.getNormal();

			Gauss3Point(v_plu3[0], v_plu3[1], v_plu3[2], vpg3);
			Gauss3Point(v_min3[0], v_min3[1], v_min3[2], vmg3);

			const auto tmp_rad = radiationFunction(vpg3, vmg3, v_plu, v_min, cen_box, arr_k_(k));
			//const auto tmp_rec = halfRecvFunction(vpg3, vmg3, v_plu, v_min, nml_plu, nml_min, center, arr_k_(k));
			//rad_mat(u, k) = (length / 2) * tmp_rad;
			//rec_mat(k, u) = (length / 2) * (conj(tmp_rad));
			rad_com_(k, u) = (length / 2) * tmp_rad;
			rec_com_(u, k) = (length / 2) * (conj(tmp_rad));
		}
	}*/

	for (int u = 0; u < com_ce; u++)
	{
		int nvp, nvm, nfp, nfm;
		VectorR3 v_plu, v_min, nml_plu, nml_min;
		VectorR3 v_plu3[3], v_min3[3], vpg3[3], vmg3[3];
		value_t length;
		ce_ptr_->getCommonEdge(u, nvp, nvm, nfp, nfm, length);
		v_plu = mesh_ptr_->getVertex(nvp);
		v_min = mesh_ptr_->getVertex(nvm);

		auto& tri_plu = mesh_ptr_->getTriangleRef(nfp);
		auto& tri_min = mesh_ptr_->getTriangleRef(nfm);
		tri_plu.getVertex(v_plu3[0], v_plu3[1], v_plu3[2]);
		tri_min.getVertex(v_min3[0], v_min3[1], v_min3[2]);
		nml_plu = tri_plu.getNormal();
		nml_min = tri_min.getNormal();

		Gauss3Point(v_plu3[0], v_plu3[1], v_plu3[2], vpg3);
		Gauss3Point(v_min3[0], v_min3[1], v_min3[2], vmg3);
		for (int k = 0; k < knum; k++)
		{
			const auto tmp_rad = radiationFunction(vpg3, vmg3, v_plu, v_min, cen_box, arr_k_(k));
			//const auto tmp_rec = halfRecvFunction(vpg3, vmg3, v_plu, v_min, nml_plu, nml_min, center, arr_k_(k));
			//rad_mat(u, k) = (length / 2) * tmp_rad;
			//rec_mat(k, u) = (length / 2) * (conj(tmp_rad));
			rad_com_(k, u) = (length / 2) * tmp_rad;
			rec_com_(u, k) = (length / 2) * (conj(tmp_rad));
		}
		dp.push_back(dipoleFunction(v_plu3, v_min3, length));
	}

	Qc3_mat rad_temp_, rec_temp_;
	rad_temp_.set_size(knum, svd_k);
	rec_temp_.set_size(svd_k, knum);
	dp_.set_size(9, svd_k);
	for (int id = 0; id < 9; id++)
	{
		for (int s_k = 0; s_k < svd_k; s_k++)
		{
			CVector I_temp = I_svd_eigen[id].col(s_k);
			for (int k = 0; k < knum; k++)
			{
				VectorC3 rad_val, rec_val;
				for (int u = 0; u < com_ce; u++)
				{
					rad_val += I_temp(u)*rad_com_(k, u);
					rec_val += I_temp(u)*rec_com_(u, k);
				}
				rad_temp_(k, s_k) = rad_val;
				rec_temp_(s_k, k) = rec_val;
			}
			VectorC3 dp_val;
			for (int u = 0; u < com_ce; u++)
				dp_val += dp[u] * I_temp(u);
			dp_(id, s_k) = dp_val;
		}
		sed_rad_.push_back(rad_temp_);
		sed_rec_.push_back(rec_temp_);
	}

}

void SED_CM_FMM::CalculateNearZ()
{
	nearZ_.resize(Array_x*Array_y*svd_k, Array_x*Array_y*svd_k);
	
	for (int idx = 0; idx < Array_x; idx++)
	{
		for (int idy = 0; idy < Array_y; idy++)
		{
			fillnearZ_kernel(idx, idy);
			/*CMatrix z_temp;
			CMatrix im_temp,in_temp;
			int idm_sed,idn_sed;
			idm_sed = ArrayToSED(idx, idy);
			im_temp = I_svd_eigen[idm_sed];

			if (idx - 1 >= 0)
			{
				//(x-1,y)
				idn_sed = ArrayToSED(idx - 1, idy);

				in_temp = I_svd_eigen[idn_sed];
				z_temp = Z_sed_eigen.block(4 * com_ce, com_ce, com_ce, com_ce);

				nearZ_.block(svd_k*(idx*Array_y + idy), svd_k*((idx - 1)*Array_y + idy), svd_k, svd_k) 
					= im_temp.transpose()*z_temp*in_temp;

				if (idy + 1 < Array_y)//(x-1,y+1)
				{
					idn_sed = ArrayToSED(idx - 1, idy+1);

					in_temp = I_svd_eigen[idn_sed];
					z_temp = Z_sed_eigen.block(4 * com_ce, 2 * com_ce, com_ce, com_ce);

					nearZ_.block(svd_k*(idx*Array_y + idy), svd_k*((idx - 1)*Array_y + idy+1), svd_k, svd_k)
						= im_temp.transpose()*z_temp*in_temp;
				}
			}

			if (idy + 1 < Array_y)
			{
				//(x,y+1)
				idn_sed = ArrayToSED(idx, idy+1);

				in_temp = I_svd_eigen[idn_sed];
				z_temp = Z_sed_eigen.block(4 * com_ce, 5*com_ce, com_ce, com_ce);

				nearZ_.block(svd_k*(idx*Array_y + idy), svd_k*(idx*Array_y + idy+1), svd_k, svd_k)
					= im_temp.transpose()*z_temp*in_temp;
				if (idx + 1 < Array_x)//(x+1,y+1)
				{
					idn_sed = ArrayToSED(idx+1, idy + 1);

					in_temp = I_svd_eigen[idn_sed];
					z_temp = Z_sed_eigen.block(4 * com_ce, 8 * com_ce, com_ce, com_ce);

					nearZ_.block(svd_k*(idx*Array_y + idy), svd_k*((idx+1)*Array_y + idy + 1), svd_k, svd_k)
						= im_temp.transpose()*z_temp*in_temp;
				}
			}

			if (idx + 1 < Array_x)
			{
				//(x+1,y)
				idn_sed = ArrayToSED(idx+1, idy);

				in_temp = I_svd_eigen[idn_sed];
				z_temp = Z_sed_eigen.block(4 * com_ce, 7 * com_ce, com_ce, com_ce);

				nearZ_.block(svd_k*(idx*Array_y + idy), svd_k*((idx + 1)*Array_y + idy), svd_k, svd_k)
					= im_temp.transpose()*z_temp*in_temp;
				if (idy - 1 >= 0)//(x+1,y-1)
				{
					idn_sed = ArrayToSED(idx + 1, idy-1);

					in_temp = I_svd_eigen[idn_sed];
					z_temp = Z_sed_eigen.block(4 * com_ce, 6 * com_ce, com_ce, com_ce);

					nearZ_.block(svd_k*(idx*Array_y + idy), svd_k*((idx + 1)*Array_y + idy-1), svd_k, svd_k)
						= im_temp.transpose()*z_temp*in_temp;
				}
			}

			if (idy - 1 >= 0)
			{
				//(x,y-1)
				idn_sed = ArrayToSED(idx, idy-1);

				in_temp = I_svd_eigen[idn_sed];
				z_temp = Z_sed_eigen.block(4 * com_ce, 3 * com_ce, com_ce, com_ce);

				nearZ_.block(svd_k*(idx*Array_y + idy), svd_k*(idx*Array_y + idy-1), svd_k, svd_k)
					= im_temp.transpose()*z_temp*in_temp;
				if (idx - 1 >= 0)//(x-1,y-1)
				{
					idn_sed = ArrayToSED(idx-1, idy - 1);

					in_temp = I_svd_eigen[idn_sed];
					z_temp = Z_sed_eigen.block(4 * com_ce, 0, com_ce, com_ce);

					nearZ_.block(svd_k*(idx*Array_y + idy), svd_k*((idx-1)*Array_y + idy - 1), svd_k, svd_k)
						= im_temp.transpose()*z_temp*in_temp;
				}
			}*/
		}
	}
}

void SED_CM_FMM::CalculateNearZ_v2()
{
	nearZ_.resize(Array_x*Array_y*svd_k, Array_x*Array_y*svd_k);
	nearZ_0inv.resize(Array_x*Array_y*svd_k, Array_x*Array_y*svd_k);
	nearZ_1.resize(Array_x*Array_y*svd_k, Array_x*Array_y*svd_k);
	
	int near_num;
	near_num = svd_k * svd_k*(16 + 6 * (2 * (Array_x - 2) + 2 * (Array_y - 2)) + 9 * ((Array_x - 2)*(Array_y - 2)));
	//nearZ_.reserve(near_num);
	//nearZ_1.reserve(near_num- Array_x * Array_y*svd_k*svd_k);
	Loc_near_ele_.resize(near_num);
	Loc_near_ele_.reserve(near_num);
	Loc_nearZ1_ele_.resize(near_num - Array_x * Array_y*svd_k*svd_k);
	Loc_nearZ1_ele_.reserve(near_num - Array_x * Array_y*svd_k*svd_k);

	auto task_num = rad_rec_assign.size() - 1;
	for (int t = 0; t < task_num; t++)
	{
		int begin = rad_rec_assign[t];
		int end = rad_rec_assign[t + 1];
		if (begin < end)
		{
			auto task = std::bind(&SED_CM_FMM::fillnearZ_byblkarr, this, begin, end);
			pool_.submit(task);
		}
	}
	pool_.run();
	nearZ_.setFromTriplets(Loc_near_ele_.begin(), Loc_near_ele_.end());
	
	if (preconditioning_)
	{
		/*Eigen::SparseLU<SpCMat> precond;
		precond.compute(nearZ_);
		CMatrix _preI;
		_preI.setIdentity(Array_x*Array_y*svd_k, Array_x*Array_y*svd_k);
		preLU_ = precond.solve(_preI);*///ILU
		nearZ_1.setFromTriplets(Loc_nearZ1_ele_.begin(), Loc_nearZ1_ele_.end());//块对角
	}
		
}

void SED_CM_FMM::preCalculateNearZ()
{
	preNearZ_.setZero(25 * svd_k, 25 * svd_k);
	for (int idx = 0; idx < 5; idx++)
	{
		for (int idy = 0; idy < 5; idy++)
		{
			fillprenearZ_kernel(idx, idy);
		}
	}
}

void SED_CM_FMM::fillnearZ_kernelMP(int begin, int end)
{
	for (int t = begin; t < end; t++)
	{
		int _idx = rad_rec_task[t].first;
		int _idy = rad_rec_task[t].second;

		CMatrix z_temp;
		CMatrix im_temp, in_temp;
		int idm_sed, idn_sed;
		idm_sed = ArrayToSED(_idx, _idy);
		im_temp = I_svd_eigen[idm_sed];

		//(x,y)
		idn_sed = ArrayToSED(_idx, _idy);

		in_temp = I_svd_eigen[idn_sed];
		z_temp = Z_sed_eigen.block(0, 0, com_ce, com_ce);

		for (int r = 0; r < svd_k; r++)
			for (int c = 0; c < svd_k; c++)
			{
				int l_r = svd_k * (_idx*Array_y + _idy) + r;
				int l_c = svd_k * (_idx*Array_y + _idy) + c;
				Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
				//Loc_near_ele_.push_back(LocSpE(l_r, l_c, val));
				{
					MutexGuard guard(mutex_);
					nearZ_.insert(l_r, l_c) = val;
				}
			}


		if (_idx - 1 >= 0)
		{
			//(x-1,y)
			idn_sed = ArrayToSED(_idx - 1, _idy);

			in_temp = I_svd_eigen[idn_sed];
			z_temp = Z_sed_eigen.block(4 * com_ce, com_ce, com_ce, com_ce);

			//nearZ_.block(svd_k*(_idx*Array_y + _idy), svd_k*((_idx - 1)*Array_y + _idy), svd_k, svd_k) = im_temp.transpose()*z_temp*in_temp;
			for (int r = 0; r < svd_k; r++)
				for (int c = 0; c < svd_k; c++)
				{
					int l_r = svd_k * (_idx*Array_y + _idy) + r;
					int l_c = svd_k * ((_idx - 1)*Array_y + _idy) + c;
					Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
					//Loc_near_ele_.push_back(LocSpE(l_r, l_c, val));
					{
						MutexGuard guard(mutex_);
						nearZ_.insert(l_r, l_c) = val;
						nearZ_1.insert(l_r, l_c) = val;
					}
				}

			if (_idy + 1 < Array_y)//(x-1,y+1)
			{
				idn_sed = ArrayToSED(_idx - 1, _idy + 1);

				in_temp = I_svd_eigen[idn_sed];
				z_temp = Z_sed_eigen.block(4 * com_ce, 2 * com_ce, com_ce, com_ce);

				//nearZ_.block(svd_k*(_idx*Array_y + _idy), svd_k*((_idx - 1)*Array_y + _idy + 1), svd_k, svd_k) = im_temp.transpose()*z_temp*in_temp;
				for (int r = 0; r < svd_k; r++)
					for (int c = 0; c < svd_k; c++)
					{
						int l_r = svd_k * (_idx*Array_y + _idy) + r;
						int l_c = svd_k * ((_idx - 1)*Array_y + _idy + 1) + c;
						Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
						//Loc_near_ele_.push_back(LocSpE(l_r, l_c, val));
						{
							MutexGuard guard(mutex_);
							nearZ_.insert(l_r, l_c) = val;
							nearZ_1.insert(l_r, l_c) = val;
						}
					}
			}
		}

		if (_idy + 1 < Array_y)
		{
			//(x,y+1)
			idn_sed = ArrayToSED(_idx, _idy + 1);

			in_temp = I_svd_eigen[idn_sed];
			z_temp = Z_sed_eigen.block(4 * com_ce, 5 * com_ce, com_ce, com_ce);

			//nearZ_.block(svd_k*(_idx*Array_y + _idy), svd_k*(_idx*Array_y + _idy + 1), svd_k, svd_k)= im_temp.transpose()*z_temp*in_temp;
			for (int r = 0; r < svd_k; r++)
				for (int c = 0; c < svd_k; c++)
				{
					int l_r = svd_k * (_idx*Array_y + _idy) + r;
					int l_c = svd_k * (_idx*Array_y + _idy + 1) + c;
					Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
					//Loc_near_ele_.push_back(LocSpE(l_r, l_c, val));
					{
						MutexGuard guard(mutex_);
						nearZ_.insert(l_r, l_c) = val;
						nearZ_1.insert(l_r, l_c) = val;
					}
				}

			if (_idx + 1 < Array_x)//(x+1,y+1)
			{
				idn_sed = ArrayToSED(_idx + 1, _idy + 1);

				in_temp = I_svd_eigen[idn_sed];
				z_temp = Z_sed_eigen.block(4 * com_ce, 8 * com_ce, com_ce, com_ce);

				//nearZ_.block(svd_k*(_idx*Array_y + _idy), svd_k*((_idx + 1)*Array_y + _idy + 1), svd_k, svd_k) = im_temp.transpose()*z_temp*in_temp;
				for (int r = 0; r < svd_k; r++)
					for (int c = 0; c < svd_k; c++)
					{
						int l_r = svd_k * (_idx*Array_y + _idy) + r;
						int l_c = svd_k * ((_idx + 1)*Array_y + _idy + 1) + c;
						Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
						//Loc_near_ele_.push_back(LocSpE(l_r, l_c, val));
						{
							MutexGuard guard(mutex_);
							nearZ_.insert(l_r, l_c) = val;
							nearZ_1.insert(l_r, l_c) = val;
						}
					}
			}
		}

		if (_idx + 1 < Array_x)
		{
			//(x+1,y)
			idn_sed = ArrayToSED(_idx + 1, _idy);

			in_temp = I_svd_eigen[idn_sed];
			z_temp = Z_sed_eigen.block(4 * com_ce, 7 * com_ce, com_ce, com_ce);

			//nearZ_.block(svd_k*(_idx*Array_y + _idy), svd_k*((_idx + 1)*Array_y + _idy), svd_k, svd_k)= im_temp.transpose()*z_temp*in_temp;
			for (int r = 0; r < svd_k; r++)
				for (int c = 0; c < svd_k; c++)
				{
					int l_r = svd_k * (_idx*Array_y + _idy) + r;
					int l_c = svd_k * ((_idx + 1)*Array_y + _idy) + c;
					Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
					//Loc_near_ele_.push_back(LocSpE(l_r, l_c, val));
					{
						MutexGuard guard(mutex_);
						nearZ_.insert(l_r, l_c) = val;
						nearZ_1.insert(l_r, l_c) = val;
					}
				}

			if (_idy - 1 >= 0)//(x+1,y-1)
			{
				idn_sed = ArrayToSED(_idx + 1, _idy - 1);

				in_temp = I_svd_eigen[idn_sed];
				z_temp = Z_sed_eigen.block(4 * com_ce, 6 * com_ce, com_ce, com_ce);

				//nearZ_.block(svd_k*(_idx*Array_y + _idy), svd_k*((_idx + 1)*Array_y + _idy - 1), svd_k, svd_k) = im_temp.transpose()*z_temp*in_temp;
				for (int r = 0; r < svd_k; r++)
					for (int c = 0; c < svd_k; c++)
					{
						int l_r = svd_k * (_idx*Array_y + _idy) + r;
						int l_c = svd_k * ((_idx + 1)*Array_y + _idy - 1) + c;
						Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
						//Loc_near_ele_.push_back(LocSpE(l_r, l_c, val));
						{
							MutexGuard guard(mutex_);
							nearZ_.insert(l_r, l_c) = val;
							nearZ_1.insert(l_r, l_c) = val;
						}
					}
			}
		}

		if (_idy - 1 >= 0)
		{
			//(x,y-1)
			idn_sed = ArrayToSED(_idx, _idy - 1);

			in_temp = I_svd_eigen[idn_sed];
			z_temp = Z_sed_eigen.block(4 * com_ce, 3 * com_ce, com_ce, com_ce);

			//nearZ_.block(svd_k*(_idx*Array_y + _idy), svd_k*(_idx*Array_y + _idy - 1), svd_k, svd_k)= im_temp.transpose()*z_temp*in_temp;
			for (int r = 0; r < svd_k; r++)
				for (int c = 0; c < svd_k; c++)
				{
					int l_r = svd_k * (_idx*Array_y + _idy) + r;
					int l_c = svd_k * (_idx*Array_y + _idy - 1) + c;
					Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
					//Loc_near_ele_.push_back(LocSpE(l_r, l_c, val));
					{
						MutexGuard guard(mutex_);
						nearZ_.insert(l_r, l_c) = val;
						nearZ_1.insert(l_r, l_c) = val;
					}
				}

			if (_idx - 1 >= 0)//(x-1,y-1)
			{
				idn_sed = ArrayToSED(_idx - 1, _idy - 1);

				in_temp = I_svd_eigen[idn_sed];
				z_temp = Z_sed_eigen.block(4 * com_ce, 0, com_ce, com_ce);

				//nearZ_.block(svd_k*(_idx*Array_y + _idy), svd_k*((_idx - 1)*Array_y + _idy - 1), svd_k, svd_k)= im_temp.transpose()*z_temp*in_temp;
				for (int r = 0; r < svd_k; r++)
					for (int c = 0; c < svd_k; c++)
					{
						int l_r = svd_k * (_idx*Array_y + _idy) + r;
						int l_c = svd_k * ((_idx - 1)*Array_y + _idy - 1) + c;
						Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
						//Loc_near_ele_.push_back(LocSpE(l_r, l_c, val));
						{
							MutexGuard guard(mutex_);
							nearZ_.insert(l_r, l_c) = val;
							nearZ_1.insert(l_r, l_c) = val;
						}
					}
			}
		}

	}
	
}

void SED_CM_FMM::fastfillnearZ_kernelMP(int begin, int end)
{
	for (int t = begin; t < end; t++)
	{
		int _idx = rad_rec_task[t].first;
		int _idy = rad_rec_task[t].second;

		int idm_near = ArrayToNear(_idx, _idy);
		//Qcout << _idx << '\t' << _idy << '\t' << idm_near << std::endl;
		int idn_near;
		int offset_m = idm_near * svd_k;
		int offset_n;
		//(x,y)
		for (int r = 0; r < svd_k; r++)
			for (int c = 0; c < svd_k; c++)
			{
				int l_r = svd_k * (_idx*Array_y + _idy) + r;
				int l_c = svd_k * (_idx*Array_y + _idy) + c;
				//nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_m + c);
				{
					MutexGuard guard(mutex_);
					nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_m + c);
				}
			}


		if (_idx - 1 >= 0)
		{
			//(x-1,y)
			//idn_near = ArrayToNear(_idx - 1, _idy);
			idn_near = idm_near - 5;
			offset_n = idn_near * svd_k;
			for (int r = 0; r < svd_k; r++)
				for (int c = 0; c < svd_k; c++)
				{
					int l_r = svd_k * (_idx*Array_y + _idy) + r;
					int l_c = svd_k * ((_idx - 1)*Array_y + _idy) + c;
					/*nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
					nearZ_1.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);*/
					{
						MutexGuard guard(mutex_);
						nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
						nearZ_1.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
					}
				}

			if (_idy + 1 < Array_y)//(x-1,y+1)
			{
				//idn_near = ArrayToNear(_idx - 1, _idy + 1);
				idn_near = idm_near - 4;
				offset_n = idn_near * svd_k;

				for (int r = 0; r < svd_k; r++)
					for (int c = 0; c < svd_k; c++)
					{
						int l_r = svd_k * (_idx*Array_y + _idy) + r;
						int l_c = svd_k * ((_idx - 1)*Array_y + _idy + 1) + c;
						/*nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
						nearZ_1.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);*/
						{
							MutexGuard guard(mutex_);
							nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
							nearZ_1.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
						}
					}
			}
		}

		if (_idy + 1 < Array_y)
		{
			//(x,y+1)
			//idn_near = ArrayToNear(_idx, _idy + 1);
			idn_near = idm_near + 1;
			offset_n = idn_near * svd_k;

			for (int r = 0; r < svd_k; r++)
				for (int c = 0; c < svd_k; c++)
				{
					int l_r = svd_k * (_idx*Array_y + _idy) + r;
					int l_c = svd_k * (_idx*Array_y + _idy + 1) + c;
					/*nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
					nearZ_1.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);*/
					{
						MutexGuard guard(mutex_);
						nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
						nearZ_1.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
					}
				}

			if (_idx + 1 < Array_x)//(x+1,y+1)
			{
				//idn_near = ArrayToNear(_idx + 1, _idy + 1);
				idn_near = idm_near + 6;
				offset_n = idn_near * svd_k;
		
				for (int r = 0; r < svd_k; r++)
					for (int c = 0; c < svd_k; c++)
					{
						int l_r = svd_k * (_idx*Array_y + _idy) + r;
						int l_c = svd_k * ((_idx + 1)*Array_y + _idy + 1) + c;
						/*nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
						nearZ_1.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);*/
						{
							MutexGuard guard(mutex_);
							nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
							nearZ_1.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
						}
					}
			}
		}

		if (_idx + 1 < Array_x)
		{
			//(x+1,y)
			//idn_near = ArrayToNear(_idx + 1, _idy);
			idn_near = idm_near + 5;
			offset_n = idn_near * svd_k;

			for (int r = 0; r < svd_k; r++)
				for (int c = 0; c < svd_k; c++)
				{
					int l_r = svd_k * (_idx*Array_y + _idy) + r;
					int l_c = svd_k * ((_idx + 1)*Array_y + _idy) + c;
					/*nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
					nearZ_1.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);*/
					{
						MutexGuard guard(mutex_);
						nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
						nearZ_1.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
					}
				}

			if (_idy - 1 >= 0)//(x+1,y-1)
			{
				//idn_near = ArrayToNear(_idx + 1, _idy - 1);
				idn_near = idm_near + 4;
				offset_n = idn_near * svd_k;

				for (int r = 0; r < svd_k; r++)
					for (int c = 0; c < svd_k; c++)
					{
						int l_r = svd_k * (_idx*Array_y + _idy) + r;
						int l_c = svd_k * ((_idx + 1)*Array_y + _idy - 1) + c;
						/*nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
						nearZ_1.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);*/
						{
							MutexGuard guard(mutex_);
							nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
							nearZ_1.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
						}
					}
			}
		}

		if (_idy - 1 >= 0)
		{
			//(x,y-1)
			//idn_near = ArrayToNear(_idx, _idy - 1);
			idn_near = idm_near - 1;
			offset_n = idn_near * svd_k;

			for (int r = 0; r < svd_k; r++)
				for (int c = 0; c < svd_k; c++)
				{
					int l_r = svd_k * (_idx*Array_y + _idy) + r;
					int l_c = svd_k * (_idx*Array_y + _idy - 1) + c;
					/*nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
					nearZ_1.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);*/
					{
						MutexGuard guard(mutex_);
						nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
						nearZ_1.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
					}
				}

			if (_idx - 1 >= 0)//(x-1,y-1)
			{
				//idn_near = ArrayToNear(_idx - 1, _idy - 1);
				idn_near = idm_near - 6;
				offset_n = idn_near * svd_k;
				
				for (int r = 0; r < svd_k; r++)
					for (int c = 0; c < svd_k; c++)
					{
						int l_r = svd_k * (_idx*Array_y + _idy) + r;
						int l_c = svd_k * ((_idx - 1)*Array_y + _idy - 1) + c;
						/*nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
						nearZ_1.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);*/
						{
							MutexGuard guard(mutex_);
							nearZ_.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
							nearZ_1.insert(l_r, l_c) = preNearZ_(offset_m + r, offset_n + c);
						}
					}
			}
		}
	}

}

void SED_CM_FMM::fillnearZ_byblkarr(int begin, int end)
{
	for (int t = begin; t < end; t++)
	{
		int idx = rad_rec_task[t].first;
		int idy = rad_rec_task[t].second;
		int id = idx * Array_y + idy;
		int offset_near_m = svd_k * BlockArray[id]._id_near;
		int cur = 0;
		int pre = 0;
		for (int r = 0; r < svd_k; r++)
		{
			for (int c = 0; c < svd_k; c++)
			{
				int l_r = svd_k * (idx*Array_y + idy) + r;
				int l_c = svd_k * (idx*Array_y + idy) + c;
				Complex val = preNearZ_(offset_near_m + r, offset_near_m + c);
				Loc_near_ele_[BlockArray[id].offset + cur] = LocSpE(l_r, l_c, val);
				cur++;
			}
		}

		for (auto &pblk : BlockArray[id].near)
		{
			int offset_near_n = svd_k * (BlockArray[id]._id_near + 5 * (pblk->_idx - idx) + pblk->_idy - idy);
			for (int r = 0; r < svd_k; r++)
			{
				for (int c = 0; c < svd_k; c++)
				{
					int l_r = svd_k * (idx*Array_y + idy) + r;
					int l_c = svd_k * (pblk->_idx*Array_y + pblk->_idy) + c;
					Complex val = preNearZ_(offset_near_m + r, offset_near_n + c);
					Loc_near_ele_[BlockArray[id].offset + cur] = LocSpE(l_r, l_c, val);
					cur++;
					Loc_nearZ1_ele_[BlockArray[id].preoffset + pre] = LocSpE(l_r, l_c, val);
					pre++;
				}
			}
		}
	}
}

void SED_CM_FMM::fillnearZ_kernel(int &_idx, int &_idy)
{
	CMatrix z_temp;
	CMatrix im_temp, in_temp;
	int idm_sed, idn_sed;
	idm_sed = ArrayToSED(_idx, _idy);
	im_temp = I_svd_eigen[idm_sed];
	
	//(x,y)
	idn_sed = ArrayToSED(_idx, _idy);

	in_temp = I_svd_eigen[idn_sed];
	z_temp = Z_sed_eigen.block(4 * com_ce, 4 * com_ce, com_ce, com_ce);

	for (int r = 0; r < svd_k; r++)
		for (int c = 0; c < svd_k; c++)
		{
			int l_r = svd_k * (_idx*Array_y + _idy) + r;
			int l_c = svd_k * (_idx*Array_y + _idy) + c;
			Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
			nearZ_.insert(l_r, l_c) = val;
		}


	if (_idx - 1 >= 0)
	{
		//(x-1,y)
		idn_sed = ArrayToSED(_idx - 1, _idy);

		in_temp = I_svd_eigen[idn_sed];
		z_temp = Z_sed_eigen.block(4 * com_ce, com_ce, com_ce, com_ce);

		//nearZ_.block(svd_k*(_idx*Array_y + _idy), svd_k*((_idx - 1)*Array_y + _idy), svd_k, svd_k) = im_temp.transpose()*z_temp*in_temp;
		for (int r = 0; r < svd_k; r++)
			for (int c = 0; c < svd_k; c++)
			{
				int l_r = svd_k * (_idx*Array_y + _idy) + r;
				int l_c = svd_k * ((_idx - 1)*Array_y + _idy) + c;
				Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
				nearZ_.insert(l_r, l_c) = val;
			}	

		if (_idy + 1 < Array_y)//(x-1,y+1)
		{
			idn_sed = ArrayToSED(_idx - 1, _idy + 1);

			in_temp = I_svd_eigen[idn_sed];
			z_temp = Z_sed_eigen.block(4 * com_ce, 2 * com_ce, com_ce, com_ce);

			//nearZ_.block(svd_k*(_idx*Array_y + _idy), svd_k*((_idx - 1)*Array_y + _idy + 1), svd_k, svd_k) = im_temp.transpose()*z_temp*in_temp;
			for (int r = 0; r < svd_k; r++)
				for (int c = 0; c < svd_k; c++)
				{
					int l_r = svd_k * (_idx*Array_y + _idy) + r;
					int l_c = svd_k * ((_idx - 1)*Array_y + _idy + 1) + c;
					Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
					nearZ_.insert(l_r, l_c) = val;
				}
		}
	}

	if (_idy + 1 < Array_y)
	{
		//(x,y+1)
		idn_sed = ArrayToSED(_idx, _idy + 1);

		in_temp = I_svd_eigen[idn_sed];
		z_temp = Z_sed_eigen.block(4 * com_ce, 5 * com_ce, com_ce, com_ce);

		//nearZ_.block(svd_k*(_idx*Array_y + _idy), svd_k*(_idx*Array_y + _idy + 1), svd_k, svd_k)= im_temp.transpose()*z_temp*in_temp;
		for (int r = 0; r < svd_k; r++)
			for (int c = 0; c < svd_k; c++)
			{
				int l_r = svd_k * (_idx*Array_y + _idy) + r;
				int l_c = svd_k * (_idx*Array_y + _idy + 1) + c;
				Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
				nearZ_.insert(l_r, l_c) = val;
			}

		if (_idx + 1 < Array_x)//(x+1,y+1)
		{
			idn_sed = ArrayToSED(_idx + 1, _idy + 1);

			in_temp = I_svd_eigen[idn_sed];
			z_temp = Z_sed_eigen.block(4 * com_ce, 8 * com_ce, com_ce, com_ce);

			//nearZ_.block(svd_k*(_idx*Array_y + _idy), svd_k*((_idx + 1)*Array_y + _idy + 1), svd_k, svd_k) = im_temp.transpose()*z_temp*in_temp;
			for (int r = 0; r < svd_k; r++)
				for (int c = 0; c < svd_k; c++)
				{
					int l_r = svd_k * (_idx*Array_y + _idy) + r;
					int l_c = svd_k * ((_idx + 1)*Array_y + _idy + 1) + c;
					Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
					nearZ_.insert(l_r, l_c) = val;
				}
		}
	}

	if (_idx + 1 < Array_x)
	{
		//(x+1,y)
		idn_sed = ArrayToSED(_idx + 1, _idy);

		in_temp = I_svd_eigen[idn_sed];
		z_temp = Z_sed_eigen.block(4 * com_ce, 7 * com_ce, com_ce, com_ce);

		//nearZ_.block(svd_k*(_idx*Array_y + _idy), svd_k*((_idx + 1)*Array_y + _idy), svd_k, svd_k)= im_temp.transpose()*z_temp*in_temp;
		for (int r = 0; r < svd_k; r++)
			for (int c = 0; c < svd_k; c++)
			{
				int l_r = svd_k * (_idx*Array_y + _idy) + r;
				int l_c = svd_k * ((_idx + 1)*Array_y + _idy) + c;
				Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
				nearZ_.insert(l_r, l_c) = val;
			}

		if (_idy - 1 >= 0)//(x+1,y-1)
		{
			idn_sed = ArrayToSED(_idx + 1, _idy - 1);

			in_temp = I_svd_eigen[idn_sed];
			z_temp = Z_sed_eigen.block(4 * com_ce, 6 * com_ce, com_ce, com_ce);

			//nearZ_.block(svd_k*(_idx*Array_y + _idy), svd_k*((_idx + 1)*Array_y + _idy - 1), svd_k, svd_k) = im_temp.transpose()*z_temp*in_temp;
			for (int r = 0; r < svd_k; r++)
				for (int c = 0; c < svd_k; c++)
				{
					int l_r = svd_k * (_idx*Array_y + _idy) + r;
					int l_c = svd_k * ((_idx + 1)*Array_y + _idy - 1) + c;
					Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
					nearZ_.insert(l_r, l_c) = val;
				}
		}
	}

	if (_idy - 1 >= 0)
	{
		//(x,y-1)
		idn_sed = ArrayToSED(_idx, _idy - 1);

		in_temp = I_svd_eigen[idn_sed];
		z_temp = Z_sed_eigen.block(4 * com_ce, 3 * com_ce, com_ce, com_ce);

		//nearZ_.block(svd_k*(_idx*Array_y + _idy), svd_k*(_idx*Array_y + _idy - 1), svd_k, svd_k)= im_temp.transpose()*z_temp*in_temp;
		for (int r = 0; r < svd_k; r++)
			for (int c = 0; c < svd_k; c++)
			{
				int l_r = svd_k * (_idx*Array_y + _idy) + r;
				int l_c = svd_k * (_idx*Array_y + _idy - 1) + c;
				Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
				nearZ_.insert(l_r, l_c) = val;
			}

		if (_idx - 1 >= 0)//(x-1,y-1)
		{
			idn_sed = ArrayToSED(_idx - 1, _idy - 1);

			in_temp = I_svd_eigen[idn_sed];
			z_temp = Z_sed_eigen.block(4 * com_ce, 0, com_ce, com_ce);

			//nearZ_.block(svd_k*(_idx*Array_y + _idy), svd_k*((_idx - 1)*Array_y + _idy - 1), svd_k, svd_k)= im_temp.transpose()*z_temp*in_temp;
			for (int r = 0; r < svd_k; r++)
				for (int c = 0; c < svd_k; c++)
				{
					int l_r = svd_k * (_idx*Array_y + _idy) + r;
					int l_c = svd_k * ((_idx - 1)*Array_y + _idy - 1) + c;
					Complex val = im_temp.transpose().row(r)*z_temp*in_temp.col(c);
					nearZ_.insert(l_r, l_c) = val;
				}
		}
	}
}

void SED_CM_FMM::fillprenearZ_kernel(int &_idx, int &_idy)
{
	CMatrix z_temp;
	CMatrix im_temp, in_temp;
	int idm_sed, idn_sed;
	idm_sed = NearArrayToSED(_idx, _idy);
	im_temp = I_svd_eigen[idm_sed].leftCols(svd_k);

	//(x,y)
	idn_sed = NearArrayToSED(_idx, _idy);

	in_temp = I_svd_eigen[idn_sed].leftCols(svd_k);
	z_temp = Z_sed_eigen.block(0, 0, com_ce, com_ce);

	preNearZ_.block(svd_k * (_idx*5 + _idy), svd_k * (_idx*5 + _idy),svd_k,svd_k)
		= im_temp.transpose()*z_temp*in_temp;

	if (_idx - 1 >= 0)
	{
		//(x-1,y)
		idn_sed = NearArrayToSED(_idx - 1, _idy);

		in_temp = I_svd_eigen[idn_sed].leftCols(svd_k);
		z_temp = Z_sed_eigen.block(4 * com_ce, com_ce, com_ce, com_ce);

		preNearZ_.block(svd_k * (_idx * 5 + _idy), svd_k * ((_idx-1) * 5 + _idy), svd_k, svd_k)
			= im_temp.transpose()*z_temp*in_temp;

		if (_idy + 1 < 5)//(x-1,y+1)
		{
			idn_sed = NearArrayToSED(_idx - 1, _idy + 1);

			in_temp = I_svd_eigen[idn_sed].leftCols(svd_k);
			z_temp = Z_sed_eigen.block(4 * com_ce, 2 * com_ce, com_ce, com_ce);

			preNearZ_.block(svd_k * (_idx * 5 + _idy), svd_k * ((_idx - 1) * 5 + _idy+1), svd_k, svd_k)
				= im_temp.transpose()*z_temp*in_temp;
		}
	}

	if (_idy + 1 < 5)
	{
		//(x,y+1)
		idn_sed = NearArrayToSED(_idx, _idy + 1);

		in_temp = I_svd_eigen[idn_sed].leftCols(svd_k);
		z_temp = Z_sed_eigen.block(4 * com_ce, 5 * com_ce, com_ce, com_ce);

		preNearZ_.block(svd_k * (_idx * 5 + _idy), svd_k * (_idx  * 5 + _idy + 1), svd_k, svd_k)
			= im_temp.transpose()*z_temp*in_temp;

		if (_idx + 1 < 5)//(x+1,y+1)
		{
			idn_sed = NearArrayToSED(_idx + 1, _idy + 1);

			in_temp = I_svd_eigen[idn_sed].leftCols(svd_k);
			z_temp = Z_sed_eigen.block(4 * com_ce, 8 * com_ce, com_ce, com_ce);

			preNearZ_.block(svd_k * (_idx * 5 + _idy), svd_k * ((_idx+1) * 5 + _idy + 1), svd_k, svd_k)
				= im_temp.transpose()*z_temp*in_temp;
		}
	}

	if (_idx + 1 < 5)
	{
		//(x+1,y)
		idn_sed = NearArrayToSED(_idx + 1, _idy);

		in_temp = I_svd_eigen[idn_sed].leftCols(svd_k);
		z_temp = Z_sed_eigen.block(4 * com_ce, 7 * com_ce, com_ce, com_ce);

		preNearZ_.block(svd_k * (_idx * 5 + _idy), svd_k * ((_idx + 1) * 5 + _idy), svd_k, svd_k)
			= im_temp.transpose()*z_temp*in_temp;

		if (_idy - 1 >= 0)//(x+1,y-1)
		{
			idn_sed = NearArrayToSED(_idx + 1, _idy - 1);

			in_temp = I_svd_eigen[idn_sed].leftCols(svd_k);
			z_temp = Z_sed_eigen.block(4 * com_ce, 6 * com_ce, com_ce, com_ce);

			preNearZ_.block(svd_k * (_idx * 5 + _idy), svd_k * ((_idx + 1) * 5 + _idy - 1), svd_k, svd_k)
				= im_temp.transpose()*z_temp*in_temp;
		}
	}

	if (_idy - 1 >= 0)
	{
		//(x,y-1)
		idn_sed = NearArrayToSED(_idx, _idy - 1);

		in_temp = I_svd_eigen[idn_sed].leftCols(svd_k);
		z_temp = Z_sed_eigen.block(4 * com_ce, 3 * com_ce, com_ce, com_ce);

		preNearZ_.block(svd_k * (_idx * 5 + _idy), svd_k * (_idx * 5 + _idy - 1), svd_k, svd_k)
			= im_temp.transpose()*z_temp*in_temp;

		if (_idx - 1 >= 0)//(x-1,y-1)
		{
			idn_sed = NearArrayToSED(_idx - 1, _idy - 1);

			in_temp = I_svd_eigen[idn_sed].leftCols(svd_k);
			z_temp = Z_sed_eigen.block(4 * com_ce, 0, com_ce, com_ce);

			preNearZ_.block(svd_k * (_idx * 5 + _idy), svd_k * ((_idx-1) * 5 + _idy - 1), svd_k, svd_k)
				= im_temp.transpose()*z_temp*in_temp;
		}
	}
}

void SED_CM_FMM::prepareSolve()
{
	const int knum = sam_phi_ * sam_theta_;
	sk_.set_size(knum, Array_x*Array_y);
	gk_.set_size(knum, Array_x*Array_y);
	dpk_.set_size(Array_x*Array_y, 1);
	coeff_.resize(knum);
	size_t cur_k = 0;
	for (int p = 0; p < sam_phi_; ++p)
		for (int t = 0; t < sam_theta_; ++t)
			coeff_(cur_k++) = gl_phi_.weight(p) * gl_theta_.weight(t) * sin(PIhalves * (gl_theta_.root(t) + 1));
	coeff_ = (k_ * k_ * Z0 / 32) * coeff_;
	coeffMatrix_ = coeff_.asDiagonal();
	assert(cur_k == knum);
}

void SED_CM_FMM::calculatePreconditioner()
{
	CMatrix z_temp = Z_sed_eigen.block(0, 0, com_ce, com_ce);
	std::vector<CMatrix> z0_inv;
	for (int id = 0; id < 9; id++)
	{
		CMatrix i_temp = I_svd_eigen[id].leftCols(svd_k);
		CMatrix z0_inv_temp = (i_temp.transpose()*z_temp*i_temp).fullPivLu().inverse();
		z0_inv.push_back(z0_inv_temp);
	}
	for (int idx = 0; idx < Array_x; idx++)
	{
		for (int idy = 0; idy < Array_y; idy++)
		{
			int id = ArrayToSED(idx, idy);
			CMatrix zinv_temp = z0_inv[id];
			int offet = (idx*Array_y + idy)*svd_k;
			for (int r = 0; r < svd_k; r++)
				for (int c = 0; c < svd_k; c++)
					nearZ_0inv.insert(offet + r, offet + c) = zinv_temp(r, c);
		}
	}
}

CVector SED_CM_FMM::matrixVectorMultiply(const int &_x, const CVector& b,const CVector& _transfer) const
{
	CVector far;
	return far;
}

CVector SED_CM_FMM::matrixVectorMultiply(const CVector& b)
{

	const int knum = sam_phi_ * sam_theta_;

	sk_.reset(VectorC3(0, 0, 0));
	for (int idx = 0; idx < Array_x; idx++)
	{
		for (int idy = 0; idy < Array_y; idy++)
		{
			int id_sed;
			id_sed = ArrayToSED(idx, idy);
			const auto & rad_mat = sed_rad_[id_sed];
			for (int r = 0; r < knum; r++)
			{
				for (int c = 0; c < svd_k; c++)
					sk_(r, idx*Array_y + idy) += rad_mat(r, c)*b((idx*Array_y + idy)*svd_k + c);
					//sk_.col(idx*Array_y + idy) = rad_mat * b.segment((idx*Array_y + idy)*svd_k, svd_k);
			}
			
		}
	}

	gk_.reset(VectorC3(0, 0, 0));
	for (int idxm = 0; idxm < Array_x; idxm++)
	{
		for (int idym = 0; idym < Array_y; idym++)
		{
			int idm = idxm * Array_y + idym;

			for (int idxn = 0; idxn < Array_x; idxn++)
			{
				for (int idyn = 0; idyn < Array_y; idyn++)
				{
					if (std::abs(idxm - idxn) > 1.5 || std::abs(idym - idyn) > 1.5)
					{
						int idn = idxn * Array_y + idyn;

						for (int k = 0; k < knum; k++)
							gk_(k, idm) += sk_(k, idn)*transfer_[k](idxm - idxn + Array_x - 1, idym - idyn + Array_y - 1);
					}
				}
			}
		}
	}

	CVector far;
	far.resize(Array_x*Array_y*svd_k);
	for (int idx = 0; idx < Array_x; idx++)
	{
		for (int idy = 0; idy < Array_y; idy++)
		{
			int id_sed,id;
			id_sed = ArrayToSED(idx, idy);
			id = idx * Array_y + idy;
			const auto & rec_mat = sed_rec_[id_sed];
			//far.segment((idx*Array_y + idy)*svd_k, svd_k) = rec_mat * coeffMatrix_* gk_.col(id);
			for (int r = 0; r < svd_k; r++)
			{
				Complex val(0, 0);
				for (int c = 0; c < knum; c++)
					val += coeff_(c)*(rec_mat(r, c)^gk_(c, id));
				far((idx*Array_y + idy)*svd_k + r) = val;
			}
		}
	}

	return nearZ_ * b + far;
}

CVector SED_CM_FMM::matrixVectorMultiply_v2(const CVector& b)
{
	CMatrix I_pre;
	I_pre.setIdentity(Array_x*Array_y*svd_k, Array_x*Array_y*svd_k);
	sk_.reset(VectorC3(0, 0, 0));
	dpk_.reset(VectorC3(0, 0, 0));
	auto task_num = rad_rec_assign.size() - 1;
	for (int t = 0; t < task_num; t++)
	{
		int begin = rad_rec_assign[t];
		int end = rad_rec_assign[t + 1];
		if (begin < end)
		{
			auto task = std::bind(&SED_CM_FMM::fillsk_kernel, this, begin, end, b);
			pool_.submit(task);
		}
	}
	pool_.run();

	gk_.reset(VectorC3(0, 0, 0));
	task_num = trans_assign.size() - 1;
	for (int t = 0; t < task_num; t++)
	{
		int begin = trans_assign[t];
		int end = trans_assign[t + 1];
		if (begin < end)
		{
			auto task = std::bind(&SED_CM_FMM::fillgk_kernel, this, begin, end);
			pool_.submit(task);
		}
	}
	pool_.run();

	farZb_.setZero(Array_x*Array_y*svd_k);
	task_num = rad_rec_assign.size() - 1;
	for (int t = 0; t < task_num; t++)
	{
		int begin = rad_rec_assign[t];
		int end = rad_rec_assign[t + 1];
		if (begin < end)
		{
			auto task = std::bind(&SED_CM_FMM::fillfarZb_kernel, this, begin, end);
			pool_.submit(task);
		}
	}
	pool_.run();

	UfarZb_.setZero(Array_x*Array_y*svd_k);
	task_num = dp_assign.size() - 1;
	for (int t = 0; t < task_num; t++)
	{
		int begin = dp_assign[t];
		int end = dp_assign[t + 1];
		if (begin < end)
		{
			auto task = std::bind(&SED_CM_FMM::fillUfarZb_kernel, this, begin, end);
			pool_.submit(task);
		}
	}
	pool_.run();

	if (preconditioning_)
	{
		return (I_pre + nearZ_01)*b + nearZ_0inv * (farZb_ + UfarZb_);//块对角
		//return preLU_ * (nearZ_ * b + farZb_ + UfarZb_);
	}
		
	return nearZ_ * b + farZb_ + UfarZb_;
}

void SED_CM_FMM::fillsk_kernel(int begin, int end,CVector &b)
{
	const int knum = sam_phi_ * sam_theta_;
	for (int t = begin; t < end; t++)
	{
		int idx = rad_rec_task[t].first;
		int idy = rad_rec_task[t].second;

		int id_sed;
		id_sed = ArrayToSED(idx, idy);
		const auto & rad_mat = sed_rad_[id_sed];
		/*for (int r = 0; r < knum; r++)
			for (int c = 0; c < svd_k; c++)
				sk_(r, idx*Array_y + idy) += rad_mat(r, c)*b((idx*Array_y + idy)*svd_k + c);*/

		for (int c = 0; c < svd_k; c++)
		{
			for (int r = 0; r < knum; r++)
				sk_(r, idx*Array_y + idy) += rad_mat(r, c)*b((idx*Array_y + idy)*svd_k + c);
			dpk_(idx*Array_y + idy, 0) += dp_(id_sed, c)*b((idx*Array_y + idy)*svd_k + c);
		}
	}
}

void SED_CM_FMM::fillgk_kernel(int begin, int end)
{
	const int knum = sam_phi_ * sam_theta_;
	for (int t = begin; t < end; t++)
	{
		int idm = trans_task[t].first;
		int idn = trans_task[t].second;
		int idxm = idm / Array_y;
		int idym = idm % Array_y;
		int idxn = idn / Array_y;
		int idyn = idn % Array_y;

		if (std::abs(idxm - idxn) > 1.5 || std::abs(idym - idyn) > 1.5)
		{
			for (int k = 0; k < knum; k++)
				gk_(k, idm) += sk_(k, idn)*transfer_[k](idxm - idxn + Array_x - 1, idym - idyn + Array_y - 1);
		}
	}
}

void SED_CM_FMM::fillfarZb_kernel(int begin, int end)
{
	const int knum = sam_phi_ * sam_theta_;
	for (int t = begin; t < end; t++)
	{
		int idx = rad_rec_task[t].first;
		int idy = rad_rec_task[t].second;

		int id_sed, id;
		id_sed = ArrayToSED(idx, idy);
		id = idx * Array_y + idy;
		const auto & rec_mat = sed_rec_[id_sed];
		for (int r = 0; r < svd_k; r++)
		{
			Complex val(0, 0);
			for (int c = 0; c < knum; c++)
				val += coeff_(c)*(rec_mat(r, c) ^ gk_(c, id));
			farZb_((idx*Array_y + idy)*svd_k + r) = val;
		}
	}
}

void SED_CM_FMM::fillUfarZb_kernel(int begin, int end)
{
	for (int t = begin; t < end; t++)
	{
		int idm = dp_task[t].first;
		int idn = dp_task[t].second;
		int idxm = idm / Array_y;
		int idym = idm % Array_y;
		int idxn = idn / Array_y;
		int idyn = idn % Array_y;

		int idm_sed = ArrayToSED(idxm, idym);
		int idn_sed = ArrayToSED(idxn, idyn);

		VectorR3 R(Dx*(idxm - idxn), Dy*(idym - idyn), 0.0f);
		VectorC3 dp_f, dp_s;

		dp_s = dpk_(idn, 0);
		for (int s_k = 0; s_k < svd_k; s_k++)
		{
			dp_f = dp_(idm_sed, s_k);
			UfarZb_((idxm*Array_y + idym)*svd_k + s_k) = dpcouplFunction(dp_f, dp_s, R);
		}

	}
}

int SED_CM_FMM::iterativeSolve_BICG()
{
	prepareSolve();
	//进行BICGSTAB迭代求解
	int iter_num = 0;
	//Qcx_vec x(unknowns_, arma::fill::zeros);
	CVector x;
	x.setZero(Array_x*Array_y*svd_k);
	//Qcx_vec r0 = V_ - matrixVectorMultiply(x);
	CVector r0 = V_red_eigen - matrixVectorMultiply(x);
	//Qcx_vec r1 = r0;
	CVector r1 = r0;
	Complex rou0(1, 0), alph(1, 0), omiga(1, 0);
	CVector v, p0, p1, s, t;
	v.setZero(x.size());
	p0.setZero(x.size());
	Complex rou1, beta;
	const auto b_normal = V_red_eigen.norm();

	auto oldState = Qcout.flags();
	Qcout << std::fixed << std::setprecision(10) << std::right;
	while (iter_num < max_iter_num_)
	{
		++iter_num;
		//rou1 = arma::cdot(r1, r0);
		rou1 = r1.dot(r0);
		beta = (rou1 / rou0) * (alph / omiga);
		p1 = r0 + beta * (p0 - omiga * v);
		v = matrixVectorMultiply(p1);
		alph = rou1 / r1.dot(v);
		s = r0 - alph * v;
		x = x + alph * p1;
		auto resdual = s.norm() / b_normal;
		Qcout << std::setw(7) << iter_num << " of " << max_iter_num_
			<< std::setw(15) << resdual << std::endl;
		if (resdual < iter_threshold)
			break;
		t = matrixVectorMultiply(s);
		//omiga = arma::cdot(t, s) / arma::cdot(t, t);
		omiga = t.dot(s) / t.dot(t);
		x = x + omiga * s;
		r0 = s - omiga * t;
		resdual = s.norm() / b_normal;
		if (resdual < iter_threshold)
			break;
		rou0 = rou1;
		p0 = p1;
	}
	I_red_eigen = x;
	Qcout.flags(oldState);
	return iter_num;
}

int SED_CM_FMM::iterativeSolve_CGN()
{
	prepareSolve();
	int iter_num = 0;
	CVector x;
	x.setZero(Array_x*Array_y*svd_k);
	CVector r0 = V_red_eigen - matrixVectorMultiply(x);
	CVector p = r0;
	//auto r0 = r.dot(r);

	auto oldState = Qcout.flags();
	Qcout << std::fixed << std::setprecision(10) << std::right;

	while (iter_num<max_iter_num_)
	{
		++iter_num;
		CVector Ap = matrixVectorMultiply_v2(p);
		Complex alpha = ((r0.transpose()*r0) / (p.transpose()*Ap))(0, 0);
		x = x + alpha * p;
		CVector r1 = r0 - alpha * Ap;

		Qcout << std::setw(7) << iter_num << " of " << max_iter_num_
			<< std::setw(15) << r1.norm() << std::endl;

		if (r1.norm() <iter_threshold)
			break;
		auto beta = ((r1.transpose()*r1) / (r0.transpose()*r0))(0, 0);
		p = r1 + beta * p;
		r0 = r1;
	}
	I_red_eigen = x;
	Qcout.flags(oldState);
	return iter_num;
}

int SED_CM_FMM::iterativeSolve_GMRES()
{
	prepareSolve();
	nearZ_01 = nearZ_0inv * nearZ_1;
	CMatrix v, h;
	CVector e, yk, x0, xk, r0;
	int m;

	v.resize(Array_x*Array_y*svd_k, max_iter_num_ + 1);
	h.setZero(max_iter_num_ + 1, max_iter_num_);
	e.setZero(max_iter_num_ + 1);
	x0.setZero(Array_x*Array_y*svd_k);
	int cur = 0;
	if (preconditioning_)
	{
		V_red_eigen = nearZ_0inv * V_red_eigen;//块对角
		//V_red_eigen = preLU_ * V_red_eigen;
	}
		

	while (true)
	{
		cur++;
		m = max_iter_num_;
		r0 = V_red_eigen - matrixVectorMultiply_v2(x0);
		auto beta = r0.norm();
		v.col(0) = r0 * (1.0f / beta);

		Qcout << std::setw(15) << cur << std::setw(15) << r0.norm() << std::endl;
		if (r0.norm() < iter_threshold)
		{
			I_red_eigen = x0;
			break;
		}

		for (int j = 0; j < max_iter_num_; j++)
		{
			CVector w = matrixVectorMultiply_v2(v.col(j));
			for (int i = 0; i < j + 1; i++)
			{
				h(i, j) = w.adjoint()*v.col(i);
				w = w - (h(i, j)*v.col(i));
			}
			//CVector v_temp = Av - v.leftCols(j)*(h.col(j).segment(0, j));
			h(j + 1, j) = w.norm();
			if (w.norm() < 1e-6)
			{
				m = j + 1;
				break;
			}
			v.col(j + 1) = w * (1.0f / h(j + 1, j));
		}
		e(0) = beta;
		yk = h.block(0, 0, m + 1, m).colPivHouseholderQr().solve(e);
		x0 = x0 + v.leftCols(m) * yk;
	}
	return cur;
}

void SED_CM_FMM::fillZbyFMM()
{
	Z_red_eigen.resize(Array_x*Array_y*svd_k, Array_x*Array_y*svd_k);
	farZ_.setZero(Array_x*Array_y*svd_k, Array_x*Array_y*svd_k);
	prepareSolve();
	int knum = sam_phi_ * sam_theta_;
	for (int idxm = 0; idxm < Array_x; idxm++)
	{
		for (int idym = 0; idym < Array_y; idym++)
		{
			auto & rec = sed_rec_[ArrayToSED(idxm, idym)];
			int idm = idxm * Array_y + idym;

			for (int idxn = 0; idxn < Array_x; idxn++)
			{
				for (int idyn = 0; idyn < Array_y; idyn++)
				{
					auto & rad = sed_rad_[ArrayToSED(idxn, idyn)];
					int idn = idxn * Array_y + idyn;
					if (std::abs(idxm - idxn) > 1.5 || std::abs(idym - idyn) > 1.5)
					{
						Qc3_mat rad_trans_;
						rad_trans_.set_size(knum, svd_k);
						rad_trans_.reset(VectorC3(0, 0, 0));
						for (int k = 0; k < knum; k++)
							for (int s_k = 0; s_k < svd_k; s_k++)
								rad_trans_(k, s_k) = rad(k, s_k)*transfer_[k](idxm - idxn + Array_x - 1, idym - idyn + Array_y - 1);
						//Qcout << "rad_trans_ is over" << std::endl;
						for (int skm = 0; skm < svd_k; skm++)
							for (int skn = 0; skn < svd_k; skn++)
								for (int k = 0; k < knum; k++)
									farZ_(idm*svd_k + skm, idn*svd_k + skn) += coeff_(k)*(rec(skm, k) ^ rad_trans_(k, skn));
						//Qcout << "farZ_ is over" << std::endl;
					}
				}
			}
		}
	}

	Z_red_eigen = nearZ_ + farZ_;
}

VectorC3 SED_CM_FMM::radiationFunction(const VectorR3* vpgs3, const VectorR3* vmgs3, const VectorR3& vp, 
	const VectorR3& vm, const VectorR3& center, const VectorR3& sam_k)
{
	VectorC3 val(0, 0, 0);
	VectorR3 rou_plu, rou_min;
	for (int p = 0; p < 3; ++p)
	{
		rou_plu = vpgs3[p] - vp;
		rou_min = vm - vmgs3[p];
		val += w3[p] * rou_plu * std::exp(J0 * k_ * (sam_k ^ (vpgs3[p] - center)));
		val += w3[p] * rou_min * std::exp(J0 * k_ * (sam_k ^ (vmgs3[p] - center)));
	}
	return val - ((sam_k ^ val) * sam_k);
}

VectorC3 SED_CM_FMM::halfRecvFunction(const VectorR3 * vpgs3, const VectorR3 * vmgs3, const VectorR3 & vp,
	const VectorR3 & vm, const VectorR3 & nml_plu, const VectorR3 & nml_min,
	const VectorR3 & center, const VectorR3 & sam_k)
{
	VectorC3 val(0, 0, 0);
	VectorR3 rou_plu, rou_min;
	for (int p = 0; p < 3; ++p)
	{
		rou_plu = vpgs3[p] - vp;
		rou_min = vm - vmgs3[p];
		val += w3[p] * (nml_plu * rou_plu) * std::exp(-J0 * k_ * (sam_k ^ (vpgs3[p] - center)));
		val += w3[p] * (nml_min * rou_min) * std::exp(-J0 * k_ * (sam_k ^ (vmgs3[p] - center)));
	}
	return sam_k * val;
}

VectorR3 SED_CM_FMM::dipoleFunction(VectorR3 *v_p3, VectorR3 *v_m3, value_t &len)
{
	VectorR3 v_cen_plu = centerv3(v_p3);
	VectorR3 v_cen_min = centerv3(v_m3);

	VectorR3 m = len * (v_cen_min - v_cen_plu);
	
	return m;
}

Complex SED_CM_FMM::dpcouplFunction(VectorC3& dp_f, VectorC3& dp_s, VectorR3& _R)
{
	value_t R_norm = _R.Norm();
	Complex C = (1.0f / (R_norm*R_norm))*(1.0f + 1.0f / (J0*k_*R_norm));
	VectorC3 M_s = (_R^dp_s)*_R / (R_norm*R_norm);

	VectorC3 E;
	E = ((M_s - dp_s)*((J0*k_ / R_norm) + C) + 2.0f*M_s*C)*exp(-J0 * k_*R_norm)*Z0 / PI4;

	return E^ dp_f;
}

void SED_CM_FMM::shiftcoordv4(int &sedid, VectorR3 *v4) const
{
	int n_x = sedid / 3;
	int n_y = sedid - n_x * 3;
	VectorR3 shiftvec(n_x*Dx, n_y*Dy, 0.0f);
	v4[0] += shiftvec;
	v4[1] += shiftvec;
	v4[2] += shiftvec;
	v4[3] += shiftvec;
}

void SED_CM_FMM::shiftcoordv1(int & sedid, VectorR3 & v) const
{
	int n_x = sedid / 3;
	int n_y = sedid - n_x * 3;
	VectorR3 shiftvec(n_x*Dx, n_y*Dy, 0.0f);
	v += shiftvec;
}

void SED_CM_FMM::shiftcoordv3(int &sedid, VectorR3 *v3) const
{
	int n_x = sedid / 3;
	int n_y = sedid - n_x * 3;
	VectorR3 shiftvec(n_x*Dx, n_y*Dy, 0.0f);
	v3[0] += shiftvec;
	v3[1] += shiftvec;
	v3[2] += shiftvec;
}

void SED_CM_FMM::shiftcoordv4(int &xid, int &yid, VectorR3 *v4) const
{
	VectorR3 shiftvec(xid*Dx, yid*Dy, 0.0f);
	v4[0] += shiftvec;
	v4[1] += shiftvec;
	v4[2] += shiftvec;
	v4[3] += shiftvec;
}

void SED_CM_FMM::shiftcoordv3(int &xid, int &yid, VectorR3 *v3) const
{
	VectorR3 shiftvec(xid*Dx, yid*Dy, 0.0f);
	v3[0] += shiftvec;
	v3[1] += shiftvec;
	v3[2] += shiftvec;
}

void SED_CM_FMM::shiftcoordv1(int & xid, int & yid, VectorR3 & v) const
{
	VectorR3 shiftvec(xid*Dx, yid*Dy, 0.0f);
	v += shiftvec;
}

void SED_CM_FMM::shiftbedgecoordv4(int &bedge, VectorR3 *v4) const
{
	if (bedge == 1)
	{
		VectorR3 shiftvec(Dx, 0.0f, 0.0f);
		v4[0] += shiftvec;
		v4[1] += shiftvec;
		v4[2] += shiftvec;
		v4[3] += shiftvec;
	}
	else
	{
		VectorR3 shiftvec(0.0f, Dy, 0.0f);
		v4[0] += shiftvec;
		v4[1] += shiftvec;
		v4[2] += shiftvec;
		v4[3] += shiftvec;
	}
}

void SED_CM_FMM::shiftbedgecoordv3(int &bedge, VectorR3 *v3) const
{
	if (bedge == 1)
	{
		VectorR3 shiftvec(Dx, 0.0f, 0.0f);
		v3[0] += shiftvec;
		v3[1] += shiftvec;
		v3[2] += shiftvec;
	}
	else
	{
		VectorR3 shiftvec(0.0f, Dy, 0.0f);
		v3[0] += shiftvec;
		v3[1] += shiftvec;
		v3[2] += shiftvec;
	}
}

void SED_CM_FMM::shiftbedgecoordv1(int & bedge, VectorR3 & v) const
{
	if (bedge == 1)
	{
		VectorR3 shiftvec(Dx, 0.0f, 0.0f);
		v += shiftvec;
	}
	else
	{
		VectorR3 shiftvec(0.0f, Dy, 0.0f);
		v += shiftvec;
	}
}

void SED_CM_FMM::shiftcoordv4(int &sedid, VectorR3 *v4, value_t &scale) const
{
	int n_x = sedid / 3;
	int n_y = sedid - n_x * 3;

	VectorR3 center_cur;
	VectorR3 shiftvec(n_x*Dx, n_y*Dy, 0.0f);

	center_cur = cen_box + shiftvec;
	v4[0] = (v4[0] - cen_box)*pow(scale, n_x) + center_cur;
	v4[1] = (v4[1] - cen_box)*pow(scale, n_x) + center_cur;
	v4[2] = (v4[2] - cen_box)*pow(scale, n_x) + center_cur;
	v4[3] = (v4[3] - cen_box)*pow(scale, n_x) + center_cur;
}

void SED_CM_FMM::shiftcoordv1(int & sedid, VectorR3 & v, value_t &scale) const
{
	int n_x = sedid / 3;
	int n_y = sedid - n_x * 3;

	VectorR3 center_cur;
	VectorR3 shiftvec(n_x*Dx, n_y*Dy, 0.0f);

	center_cur = cen_box + shiftvec;
	v = (v - cen_box)*pow(scale, n_x) + center_cur;
}

void SED_CM_FMM::shiftcoordv3(int &sedid, VectorR3 *v3, value_t &scale) const
{
	int n_x = sedid / 3;
	int n_y = sedid - n_x * 3;

	VectorR3 center_cur;
	VectorR3 shiftvec(n_x*Dx, n_y*Dy, 0.0f);

	center_cur = cen_box + shiftvec;
	v3[0] = (v3[0] - cen_box)*pow(scale, n_x) + center_cur;
	v3[1] = (v3[1] - cen_box)*pow(scale, n_x) + center_cur;
	v3[2] = (v3[2] - cen_box)*pow(scale, n_x) + center_cur;
}

void SED_CM_FMM::shiftcoordv4(int &xid, int &yid, VectorR3 *v4, value_t &scale) const
{
	VectorR3 center_cur;
	VectorR3 shiftvec(xid*Dx, yid*Dy, 0.0f);

	center_cur = cen_box + shiftvec;
	v4[0] = (v4[0] - cen_box)*pow(scale, xid) + center_cur;
	v4[1] = (v4[1] - cen_box)*pow(scale, xid) + center_cur;
	v4[2] = (v4[2] - cen_box)*pow(scale, xid) + center_cur;
	v4[3] = (v4[3] - cen_box)*pow(scale, xid) + center_cur;
}

void SED_CM_FMM::shiftcoordv3(int &xid, int &yid, VectorR3 *v3, value_t &scale) const
{
	VectorR3 center_cur;
	VectorR3 shiftvec(xid*Dx, yid*Dy, 0.0f);

	center_cur = cen_box + shiftvec;
	v3[0] = (v3[0] - cen_box)*pow(scale, xid) + center_cur;
	v3[1] = (v3[1] - cen_box)*pow(scale, xid) + center_cur;
	v3[2] = (v3[2] - cen_box)*pow(scale, xid) + center_cur;
}

void SED_CM_FMM::shiftcoordv1(int & xid, int & yid, VectorR3 & v, value_t &scale) const
{
	VectorR3 center_cur;
	VectorR3 shiftvec(xid*Dx, yid*Dy, 0.0f);

	center_cur = cen_box + shiftvec;
	v = (v - cen_box)*pow(scale, xid) + center_cur;
}

void SED_CM_FMM::scalecoordv1(int & xid, int & yid, VectorR3 & v, value_t &scale) const
{
	v = (v - cen_box)*pow(scale, xid) + cen_box;
}

void SED_CM_FMM::scalecoordv3(int & xid, int & yid, VectorR3 *v3, value_t &scale) const
{
	v3[0] = (v3[0] - cen_box)*pow(scale, xid) + cen_box;
	v3[1] = (v3[1] - cen_box)*pow(scale, xid) + cen_box;
	v3[2] = (v3[2] - cen_box)*pow(scale, xid) + cen_box;
}

void SED_CM_FMM::shiftbedgecoordv4(int &bedge, VectorR3 *v4, value_t &scale) const
{
	if (bedge == 1)
	{
		VectorR3 shiftvec(Dx, 0.0f, 0.0f);
		v4[0] += shiftvec;
		v4[1] += shiftvec;
		v4[2] += shiftvec;
		v4[3] += shiftvec;
	}
	else
	{
		VectorR3 shiftvec(0.0f, Dy, 0.0f);
		v4[0] += shiftvec;
		v4[1] += shiftvec;
		v4[2] += shiftvec;
		v4[3] += shiftvec;
	}
}

void SED_CM_FMM::shiftbedgecoordv3(int &bedge, VectorR3 *v3, value_t &scale) const
{
	if (bedge == 1)
	{
		VectorR3 shiftvec(Dx, 0.0f, 0.0f);
		v3[0] += shiftvec;
		v3[1] += shiftvec;
		v3[2] += shiftvec;
	}
	else
	{
		VectorR3 shiftvec(0.0f, Dy, 0.0f);
		v3[0] += shiftvec;
		v3[1] += shiftvec;
		v3[2] += shiftvec;
	}
}

void SED_CM_FMM::shiftbedgecoordv1(int & bedge, VectorR3 & v, value_t &scale) const
{
	if (bedge == 1)
	{
		VectorR3 shiftvec(Dx, 0.0f, 0.0f);
		v += shiftvec;
	}
	else
	{
		VectorR3 shiftvec(0.0f, Dy, 0.0f);
		v += shiftvec;
	}
}

int SED_CM_FMM::ArrayToSED(int _x, int _y) const
{
	if (_x == 0)
	{
		if (_y == 0)
		{
			//左上
			return (0);
		}
		else if (_y == (Array_y - 1))
		{
			//右上
			return (2);
		}
		else
		{
			//上边
			return (1);
		}
	}
	else if (_x == (Array_x - 1))
	{
		if (_y == 0)
		{
			//左下
			return (6);
		}
		else if (_y == (Array_y - 1))
		{
			//右下
			return (8);
		}
		else
		{
			//下边
			return (7);
		}
	}
	else
	{
		if (_y == 0)
		{
			//左边
			return (3);
		}
		else if (_y == (Array_y - 1))
		{
			//右边
			return (5);
		}
		else
		{
			//中间
			return (4);
		}
	}
}

int SED_CM_FMM::ArrayToNear(int _x, int _y) const
{
	if (Array_y > 4 && Array_x > 4)
		if (_x > 1 && _x < (Array_x - 2) && _y>1 && _y < (Array_y - 2)) return (12);

	if (Array_y > 4)
	{
		if (_x == 0 && _y > 1 && _y < ((Array_y - 2))) return (2);
		if (_x == 1 && _y > 1 && _y < ((Array_y - 2))) return (7);
		if (_x == (Array_x - 2) && _y > 1 && _y < ((Array_y - 2))) return (17);
		if (_x == (Array_x - 1) && _y > 1 && _y < ((Array_y - 2))) return (22);
	}

	if (Array_x > 4)
	{
		if (_y == 0 && _x > 1 && _x < (Array_x - 2)) return (10);
		if (_y == 1 && _x > 1 && _x < (Array_x - 2)) return (11);
		if (_y == (Array_y - 1) && _x > 1 && _x < (Array_x - 2)) return (14);
		if (_y == (Array_y - 2) && _x > 1 && _x < (Array_x - 2)) return (13);
	}

	if (_x == 0 && _y == 0) return (0);
	if (_x == 0 && _y == (Array_y - 1)) return (4);
	if (_x == 1 && _y == 0) return (5);
	if (_x == 1 && _y == (Array_y - 1)) return (9);
	if (_x == (Array_x - 2) && _y == 0) return (15);
	if (_x == (Array_x - 2) && _y == (Array_y - 1)) return (19);
	if (_x == (Array_x - 1) && _y == 0) return (20);
	if (_x == (Array_x - 1) && _y == (Array_y - 1)) return (24);
	
	if (_x == 0 && _y == 1) return (1);
	if (_x == 0 && _y == (Array_y - 2)) return(3);
	if (_x == 1 && _y == 1) return (6);
	if (_x == 1 && _y == (Array_y - 2)) return(8);
	if (_x == (Array_x - 2) && _y == 1) return (16);
	if (_x == (Array_x - 2) && _y == (Array_y - 2)) return(18);
	if (_x == (Array_x - 1) && _y == 1) return (21);
	if (_x == (Array_x - 1) && _y == (Array_y - 2)) return(23);
}

int SED_CM_FMM::NearArrayToSED(int _x, int _y) const
{
	if (_x == 0)
	{
		if (_y == 0)
		{
			//左上
			return (0);
		}
		else if (_y == (5 - 1))
		{
			//右上
			return (2);
		}
		else
		{
			//上边
			return (1);
		}
	}
	else if (_x == (5 - 1))
	{
		if (_y == 0)
		{
			//左下
			return (6);
		}
		else if (_y == (5 - 1))
		{
			//右下
			return (8);
		}
		else
		{
			//下边
			return (7);
		}
	}
	else
	{
		if (_y == 0)
		{
			//左边
			return (3);
		}
		else if (_y == (5 - 1))
		{
			//右边
			return (5);
		}
		else
		{
			//中间
			return (4);
		}
	}
}

void SED_CM_FMM::getBoundaryNum(int &_sedid, int &_b1, int &_b2) const
{
	if (_sedid == 0 || _sedid == 1 || _sedid == 3 || _sedid == 4)
	{
		_b1 = bou_ce1;
		_b2 = bou_ce2;
	}
	else if (_sedid == 2 || _sedid == 5)
	{
		_b1 = bou_ce1;
		_b2 = 0;
	}
	else if (_sedid == 6 || _sedid == 7)
	{
		_b1 = 0;
		_b2 = bou_ce2;
	}
	else
	{
		_b1 = 0;
		_b2 = 0;
	}
}