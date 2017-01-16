/*****************************
  The program requires the following files to run properly:
  a. Real location
  b. odometry information
  c. corresponding measurements (for environmental field, 0 - temperature, 1 - humid, 2 - light, 3 - others)
     FIELD names: REAL, WIRELESS, ENVIR
  d. training data: inputs, outputs, hypers
  ****************************/
#include <cstdio>
#include <cstring>
#include <cmath>
#include "nr3.h"
#include "ran.h"
#include "deviates.h"
#include "mvas_gp.h"
#include "fullgp.h"
#define _USE_MATH_DEFINES

const int STR_LEN = 256;
const int FIELD_NUM = 10;
const int TRAIN_NUM = 2024;
const int HYPER_NUM = 8;
const int STEP_NUM = 2024;
const int PARTICLE_NUM = 810;
const int ALGO_NUM = 8;

class GPBF
{
public:
	void parse(char * conf_name);
	void init();
	void readData();
	void sampleTraj(int t);
	double regulate(double v);
	VecDoub randomPose(VecDoub mp);
	VecDoub motionModel(VecDoub lpose, VecDoub lodom, VecDoub codom);
	VecDoub dMotionModel(VecDoub lpose, VecDoub lodom, VecDoub codom);
	void filter();
	void pFilter(int algo_id);
	void getSupSet();

private:
	// 0 - trunc, 1 - SoD, 2 - OPITC, 3 - IOPITC, 4 - FGP, 5 - PITC
	double alpha1, alpha2, alpha3, alpha4;
	char output_name[STR_LEN];
	int SUPPORT_NUM, SLICE_SIZE, TRUNC_NUM, case_id;
	FullGP fgp[FIELD_NUM];
	int step_num, field_num, par_num, algo_num, algo[ALGO_NUM];
	int train_num[FIELD_NUM], hyper_num[FIELD_NUM];
	int fid_list[FIELD_NUM];		//fid specifies the particular map to be used (refer to c)

	VecDoub train[FIELD_NUM][TRAIN_NUM];		//sampled training data
	VecDoub odom[STEP_NUM], real_loc[STEP_NUM], obs[STEP_NUM];
	VecDoub es_loc[STEP_NUM];
	/***********
	obs and real_loc all starts at 1, for obs, the 2nd dim indicates the measurement from each field
	odom length is step_num+1, the initial odom (0, 0, 0) stores in step_0
	************/
	VecDoub traj[PARTICLE_NUM][STEP_NUM];
	//trajectory length is step_num+1, step_0 is the initial sample pose
	char field_name[STR_LEN];
	VecDoub hyp[FIELD_NUM];
	VecDoub bel[STEP_NUM][PARTICLE_NUM];
	VecDoub bel_bar[PARTICLE_NUM];
	VecDoub pmean[PARTICLE_NUM][FIELD_NUM], pvar[PARTICLE_NUM][FIELD_NUM];
	//3rd dim is the pmean of all the particles

	/*************************
	 following variables defined for online GP (aka, SoD, OPITC, IOPITC)
	 *************************/
	bool update_time[STEP_NUM];								//for SoD only
	VecDoub sod_train[PARTICLE_NUM][FIELD_NUM][TRAIN_NUM];	//for SoD only
	int sod_train_num;										//for SoD only
	//global summary
	MatDoub gs_kuu[PARTICLE_NUM][FIELD_NUM];
    VecDoub gs_zu[PARTICLE_NUM][FIELD_NUM];
    Cholesky *chol_kuu[FIELD_NUM];
    Cholesky *chol_suu[PARTICLE_NUM][FIELD_NUM];
    //support set related
    int uset_num;
    MatDoub uset;	//each sample trajectory has a support set
    MatDoub kuu[FIELD_NUM];
    struct mvas_gp gp[FIELD_NUM];
    //for improved online PITC
    int ldata_cnt;
    MatDoub ldata[PARTICLE_NUM][FIELD_NUM], test_data;
    VecDoub lmean[PARTICLE_NUM][FIELD_NUM], lvar;		//dummy lvar

    void deadReckon();
    void initPITC();
    void onlinePITC(VecDoub &weight, int t, int flag);
    void GlobalSum(int t);
    void IOPITC_Sum(int t);
    void IOPITC_Regr(int t);
    void FGP(VecDoub &weight, int t);
	void trunc(VecDoub &weight, int t);
	void initSoD();
	void SoD(VecDoub &weight, int t);
	void batchPITC(VecDoub &weight, int t);
	void releaseMem();
	void getObs(VecDoub &weight, int t)
    {
    	for (int idx = 0; idx < par_num; idx++)
    	{
    		double expectation = 0;

	        for (int i = 0; i < par_num; i++)
	        {
	            double tv = 1;
	            for (int f = 0; f < field_num; f++)
	            {
	                double tmp = -1.0 * pow(obs[t][f] - pmean[i][f][idx], 2) / (2.0 * pvar[i][f][idx]);
	                tv *= exp(tmp) / sqrt(2.0 * M_PI * pvar[i][f][idx]);
	                // printf("%lf %lf obs=%lf tmp=%lf v=%lf\n", pmean[i][f][idx], pvar[i][f][idx], obs[t][f], tmp, tv);
	            }
	            expectation += pow(tv, 1.0 / field_num);
	        }
	        weight[idx] = expectation / par_num;
	    }
    }
};


void GPBF::parse(char * conf_name)
{
	FILE * fp_conf = fopen(conf_name, "r");
	if (fp_conf == NULL)
	{
		printf("error opening config file!\n");
		return;
	}

	// parse the config file
	char str[STR_LEN];
	while (fscanf(fp_conf, "%s", str) == 1)
	{
		if (strcmp(str, "step_num") == 0)
			fscanf(fp_conf, "%d", &step_num);
		else if (strcmp(str, "field_num") == 0)
		{
			// read the id of each field
			fscanf(fp_conf, "%s", str);
			field_num = str[0] - '0';
			printf("field number: %d\n", field_num);
			for (int i = 0; i < field_num; i++)
				fid_list[i] = str[2*(1+i)] - '0';
		}
		else if (strcmp(str, "field_name") == 0)
			fscanf(fp_conf, "%s", field_name);
		else if (strcmp(str, "particle_num") == 0)
			fscanf(fp_conf, "%d", &par_num);
		else if (strcmp(str, "algo") == 0)
		{
			fscanf(fp_conf, "%s", str);
			algo_num = str[0] - '0';
			for (int i = 0; i < algo_num; i++)
			{
				algo[i] = str[2*(i+1)]-'0';
				printf("algorithm #%d\n", algo[i]);
			}
		}
		else if (strcmp(str, "alpha") == 0)
		{
			fscanf(fp_conf, "%lf %lf %lf %lf", &alpha1, &alpha2, &alpha3, &alpha4);
			printf("alpha: %lf %lf %lf %lf\n", alpha1, alpha2, alpha3, alpha4);
		}
		else if (strcmp(str, "support_num") == 0)
			fscanf(fp_conf, "%d", &SUPPORT_NUM);
		else if (strcmp(str, "slice_size") == 0)
			fscanf(fp_conf, "%d", &SLICE_SIZE);
		else if (strcmp(str, "trunc_size") == 0)
			fscanf(fp_conf, "%d", &TRUNC_NUM);
		else if (strcmp(str, "dataset") == 0)
			fscanf(fp_conf, "%d", &case_id);
		else if (strcmp(str, "output") == 0)
			fscanf(fp_conf, "%s", output_name);
		else
		{
			printf("unexpected identifier %s!\n", str);
			return;
		}
	}
	fclose(fp_conf);
}

/********************
 * obtain the training data - input, output and hyper
 * note that the name does't have to do with field name, different trials share the same training data
 ********************/
void GPBF::init()
{
	char str[STR_LEN];
	double dummy;	//dummy variable

	srand(time(0));
	for (int i = 0; i < field_num; i++)
	{
		//read hyper
		FILE * fp_hyp = NULL;
		sprintf(str, "../dom%d/hyper_%d", case_id, fid_list[i]);
		fp_hyp = fopen(str, "r");
		if (fp_hyp == NULL)
			printf("error opening %s!\n", str);
		fscanf(fp_hyp, "%d", &hyper_num[i]);
		hyp[i].resize(hyper_num[i]);
		for (int j = 0; j < hyper_num[i]; j++)
			fscanf(fp_hyp, "%lf", &hyp[i][j]);
		fclose(fp_hyp);

		//read input and output
		FILE * fp_in = NULL;
		FILE * fp_out = NULL;
		sprintf(str, "../dom%d/output_%d", case_id, fid_list[i]);
		fp_out = fopen(str, "r");
		if (fp_out == NULL)
			printf("error opening %s!\n", str);
		//get the size of training data
		train_num[i] = 0;
		while (fscanf(fp_out, "%lf", &dummy) == 1)
			train_num[i]++;
		fclose(fp_out);

		//now get the training data, first 2 dim indicates the location, the third dim is the measurement
		sprintf(str, "../dom%d/input_%d", case_id, fid_list[i]);
		fp_in = fopen(str, "r");
		sprintf(str, "../dom%d/output_%d", case_id, fid_list[i]);
		fp_out = fopen(str, "r");
		for (int j = 0; j < train_num[i]; j++)
		{
			train[i][j].resize(3);
			fscanf(fp_in, "%lf %lf", &train[i][j][0], &train[i][j][1]);
			fscanf(fp_out, "%lf", &train[i][j][2]);
		}
		fclose(fp_in);
		fclose(fp_out);
	}
	for (int i = 0; i < field_num; i++)
		fgp[i].init(hyp[i]);
}

//obtain the real location, odometry and measurement
void GPBF::readData()
{
	FILE * fp_odom = NULL;
	FILE * fp_obs = NULL;
	FILE * fp_loc = NULL;
	char str[STR_LEN];

	sprintf(str, "../traj%d/%s_odom", case_id, field_name);
	fp_odom = fopen(str, "r");
	if (fp_odom == NULL)
		printf("error opening %s!\n", str);
	sprintf(str, "../traj%d/%s_loc", case_id, field_name);
	fp_loc = fopen(str, "r");
	if (fp_loc == NULL)
		printf("error opening %s!\n", str);


	for (int i = 0; i <= step_num; i++)		// the inital position is in step 0
	{
		real_loc[i].resize(3);
		odom[i].resize(3);

		fscanf(fp_odom, "%lf %lf %lf", &odom[i][0], &odom[i][1], &odom[i][2]);
		fscanf(fp_loc, "%lf %lf %lf", &real_loc[i][0], &real_loc[i][1], &real_loc[i][2]);
	}

	// if don't have observation file then do the prediction online (only real data have the obs file)
	sprintf(str, "../traj%d/%s_obs", case_id, field_name);
	fp_obs = fopen(str, "r");
	for (int i = 1; i <= step_num; i++)
		obs[i].resize(field_num);
	if (fp_obs == NULL)
	{
		VecDoub tmean(STEP_NUM), tvar(STEP_NUM), ttest[STEP_NUM];
		for (int i = 1; i <= step_num; i++)
		{
			ttest[i-1] = real_loc[i];
		}

		for (int f = 0; f < field_num; f++)
		{
			fgp[f].predict(train[f], train_num[f], ttest, tmean, tvar, step_num);
			for (int i = 1; i <= step_num; i++)
			{
				obs[i][f] = tmean[i-1];
			}
		}
	}
	else
	{
		double dummy_arr[3];
		for (int i = 1; i <= step_num; i++)
		{
			fscanf(fp_obs, "%lf %lf %lf", &dummy_arr[0], &dummy_arr[1], &dummy_arr[2]);
			for (int j = 0; j < field_num; j++)
			{
				obs[i][j] = dummy_arr[fid_list[j]];
			}
		}

		fclose(fp_obs);
	}


	fclose(fp_odom);
	fclose(fp_loc);
}

//deterministic version
VecDoub GPBF::dMotionModel(VecDoub lpose, VecDoub lodom, VecDoub codom)
{
	VecDoub cpose(3);

	double rot1 = regulate(atan2(codom[1] - lodom[1], codom[0] - lodom[0]) - lodom[2]);
	double trans = sqrt(pow(codom[0] - lodom[0], 2) + pow(codom[1] - lodom[1], 2));
	double rot2 = regulate(codom[2] - lodom[2] - rot1);

	cpose[0] = lpose[0] + trans * cos(lpose[2] + rot1);
	cpose[1] = lpose[1] + trans * sin(lpose[2] + rot1);
	cpose[2] = regulate(lpose[2] + rot1 + rot2);

	return cpose;
}

bool judge(VecDoub a, MatDoub mat, int num)
{
	for (int i = 0; i < num; i++)
	{
		double dist = pow(a[0] - mat[i][0], 2) + pow(a[1] - mat[i][1], 2);
		dist = sqrt(dist);
		if (dist < 0.5)
			return false;
	}
	return true;
}
void GPBF::getSupSet()
{
	es_loc[0] = real_loc[0];
	for (int i = 1; i <= step_num; i++)
	{
		update_time[i] = false;
		es_loc[i] = dMotionModel(es_loc[i-1], odom[i-1], odom[i]);
	}
	//initialize the support set for each trajectory
	uset_num = 0;
    uset.resize(SUPPORT_NUM, 3);

    while (uset_num < SUPPORT_NUM)
    {
        int tmp = rand() % step_num + 1;
        if (!judge(es_loc[tmp], uset, uset_num))
            continue;

        uset[uset_num][0] = es_loc[tmp][0];
        uset[uset_num][1] = es_loc[tmp][1];
        update_time[tmp] = true;
        uset_num++;
    }
}

//motion model
VecDoub GPBF::motionModel(VecDoub lpose, VecDoub lodom, VecDoub codom)
{
	VecDoub cpose(3);

	double rot1 = regulate(atan2(codom[1] - lodom[1], codom[0] - lodom[0]) - lodom[2]);
	double trans = sqrt(pow(codom[0] - lodom[0], 2) + pow(codom[1] - lodom[1], 2));
	double rot2 = regulate(codom[2] - lodom[2] - rot1);

	Normaldev gauss1(0, sqrt(alpha1 * rot1 * rot1 + alpha2 * trans * trans), rand());
	double hat_rot1 = rot1 - gauss1.dev();
	Normaldev gauss2(0, sqrt(alpha3 * trans * trans + alpha4 * (rot1 * rot1 + rot2 * rot2)), rand());
	double hat_trans = trans - gauss2.dev();
	Normaldev gauss3(0, sqrt(alpha1 * rot2 * rot2 + alpha2 * trans * trans), rand());
	double hat_rot2 = rot2 - gauss3.dev();

	cpose[0] = lpose[0] + hat_trans * cos(lpose[2] + hat_rot1);
	cpose[1] = lpose[1] + hat_trans * sin(lpose[2] + hat_rot1);
	cpose[2] = regulate(lpose[2] + hat_rot1 + hat_rot2);

	return cpose;
}

double GPBF::regulate(double v)
{
	while (fabs(v) > M_PI)
	{
		if (v > 0)
			v -= 2 * M_PI;
		else
			v += 2 * M_PI;
	}

	return v;
}

VecDoub GPBF::randomPose(VecDoub mp)
{
	Normaldev gauss_x(mp[0], 2, rand());
	Normaldev gauss_y(mp[1], 2, rand());
	Normaldev gauss_rad(mp[2], 0.1, rand());
	VecDoub r(3);

	r[0] = gauss_x.dev();
	r[1] = gauss_y.dev();
	r[2] = regulate(gauss_rad.dev());

	return r;
}

//generate sample trajectory
void GPBF::sampleTraj(int t)
{
	for (int i = 0; i < par_num; i++)
	{
		for (int j = t; j < t + SLICE_SIZE && j <= step_num; j++)
		{
			if (j == t)
				traj[i][j] = motionModel(bel[j-1][i], odom[j-1], odom[j]);
			else
				traj[i][j] = motionModel(traj[i][j-1], odom[j-1], odom[j]);
		}
	}
}

/*********************
 * truncation method, use the most recent data of size TRUNC_NUM
 *********************/
void GPBF::trunc(VecDoub &weight, int t)
{
	int start_idx = t > TRUNC_NUM ? t - TRUNC_NUM : 1;
	VecDoub *ttrain = new VecDoub[TRUNC_NUM];

	for (int i = 0; i < TRUNC_NUM; i++)
		ttrain[i].resize(3);
	for (int s = 0; s < par_num; s++)
	{
		for (int f = 0; f < field_num; f++)
		{
			for (int i = start_idx; i < t; i++)
			{
				ttrain[i-start_idx] = traj[s][i];
				ttrain[i-start_idx][2] = obs[i][f];
			}
			fgp[f].predict(ttrain, t - start_idx, bel_bar, pmean[s][f], pvar[s][f], par_num);
		}
	}
	getObs(weight, t);
	if ((t - 1) % SLICE_SIZE == 0)
		sampleTraj(t);
	delete[] ttrain;
}

/*********************
 * full GP, use all the past data (start from 1)
 *********************/
void GPBF::FGP(VecDoub &weight, int t)
{
	VecDoub ttrain[STEP_NUM];

	for (int i = 0; i < step_num; i++)
		ttrain[i].resize(3);
	for (int s = 0; s < par_num; s++)
	{
		for (int f = 0; f < field_num; f++)
		{
			for (int i = 1; i < t; i++)
			{
				ttrain[i-1][0] = traj[s][i][0];
				ttrain[i-1][1] = traj[s][i][1];
				ttrain[i-1][2] = obs[i][f];
			}
			fgp[f].predict(ttrain, t - 1, bel_bar, pmean[s][f], pvar[s][f], par_num);
		}
	}
	getObs(weight, t);
	if ((t - 1) % SLICE_SIZE == 0)
		sampleTraj(t);
}

/*********************
 * subset of data (SoD), use a subset of data for regression
 *********************/
void GPBF::initSoD()
{
	sod_train_num = 0;
	for (int i = 0; i < par_num; i++)
		for (int f = 0; f < field_num; f++)
			for (int j = 0; j < SUPPORT_NUM; j++)
				sod_train[i][f][j].resize(3);
}
void GPBF::SoD(VecDoub &weight, int t)
{
	for (int s = 0; s < par_num; s++)
	{
		for (int f = 0; f < field_num; f++)
		{
			fgp[f].predict(sod_train[s][f], sod_train_num, bel_bar, pmean[s][f], pvar[s][f], par_num);
		}
	}
	getObs(weight, t);
	if ((t - 1) % SLICE_SIZE == 0)
		sampleTraj(t);
	if (update_time[t])
	{
		for (int i = 0; i < par_num; i++)
		{
			for (int f = 0; f < field_num; f++)
			{
				sod_train[i][f][sod_train_num][0] = traj[i][t][0];
				sod_train[i][f][sod_train_num][1] = traj[i][t][1];
				sod_train[i][f][sod_train_num][2] = obs[t][f];
			}
		}
		sod_train_num++;
	}
}

/*********************
 * online PITC and improved online PITC
 *********************/
void GPBF::initPITC()
{
    //init gp
    for (int i = 0; i < field_num; i++)
    {
    	gp[i].init(hyp[i], 2);    //2 is the input dimension
    	for (int j = 0; j < hyp[i].size(); j++)
    		printf("%lf ", hyp[i][j]);
    	printf("\n");
    }


    //initialize global summary
    for (int f = 0; f < field_num; f++)
    {
    	gp[f].cov->se_ard(uset, uset_num, kuu[f]);
        for (int i = 0; i < par_num; i++)
        {
            gs_kuu[i][f].resize(kuu[f].nrows(), kuu[f].ncols());
            gs_kuu[i][f] = kuu[f];
            gs_zu[i][f].assign(uset_num, 0);
        }
    }

    for (int k = 0; k < field_num; k++)
    {
    	chol_kuu[k] = gp[k].chol_cov(kuu[k]);
    	for (int i = 0; i < par_num; i++)
    		chol_suu[i][k] = gp[k].chol_cov(gs_kuu[i][k]);
    }
    printf("obtained support set successfully\n");

    //init lmean, lvar, pmean and pvar and the backup
    for (int i = 0; i < par_num; i++)
    {
        for (int k = 0; k < field_num; k++)
        {
        	lmean[i][k].resize(SLICE_SIZE);
            for (int j = 0; j < par_num; j++)
            {
                pmean[i][k][j] = hyp[k][hyper_num[k]-1];
                pvar[i][k][j] = hyp[k][0] * hyp[k][0] + hyp[k][hyper_num[k]-2] * hyp[k][hyper_num[k]-2];
            }
            for (int j = 0; j < SLICE_SIZE; j++)
            	lmean[i][k][j] = hyp[k][hyper_num[k]-1];
        }
    }

    test_data.resize(par_num, 3);
    lvar.resize(SLICE_SIZE);

    //init the data for storing the data within each slice
	ldata_cnt = 1;
	for (int i = 0; i < par_num; i++)
		for (int j = 0; j < field_num; j++)
		{
			ldata[i][j].resize(SLICE_SIZE, 3);
			VecDoub tpose(3), ppose(3);
            for (int r = 1; r <= SLICE_SIZE && r <= step_num; r++)
            {
            	if (r == 1)
					tpose = motionModel(bel[0][i], odom[r-1], odom[r]);
				else
					tpose = motionModel(ppose, odom[r-1], odom[r]);
            	ldata[i][j][r-1][0] = tpose[0];
	            ldata[i][j][r-1][1] = tpose[1];
	            ldata[i][j][r-1][2] = obs[r][j];
	            ppose = tpose;
            }
		}



    printf("init pitc done\n");
}
void GPBF::GlobalSum(int t)
{
    MatDoub ls_kuu;
    VecDoub ls_zu;
    VecDoub traj_pitc[PARTICLE_NUM][SLICE_SIZE];

    MatDoub D;
    D.resize(SLICE_SIZE, 3);

    for (int i = 0; i < par_num; i++)
	{
		for (int j = t - SLICE_SIZE + 1; j <= t; j++)
		{
			if (j == t - SLICE_SIZE + 1)
				traj_pitc[i][j-t+SLICE_SIZE-1] = motionModel(bel[j-1][i], odom[j-1], odom[j]);
			else
				traj_pitc[i][j-t+SLICE_SIZE-1] = motionModel(traj_pitc[i][j-t+SLICE_SIZE-2], odom[j-1], odom[j]);
		}
	}

    for (int k = 0; k < field_num; k++)
    {
        for (int i = 0; i < par_num; i++)
        {
            for (int j = 0; j < SLICE_SIZE; j++)      // summarize data from t - SLICE_SIZE + 1 to t
            {
                //init train data
                D[j][0] = traj_pitc[i][j][0];
                D[j][1] = traj_pitc[i][j][1];
                D[j][2] = obs[j+t-SLICE_SIZE+1][k];
            }

            //local summary in block i
            Cholesky *chol_sdd = gp[k].chol_pcov(D, SLICE_SIZE, uset, uset_num, chol_kuu[k]);
            gp[k].pitc_prep(D, SLICE_SIZE, chol_sdd, uset, uset_num, chol_kuu[k], ls_zu, ls_kuu);

            //global summary
            for(Int r = 0; r < uset_num; r++)
            {
                gs_zu[i][k][r] += ls_zu[r];
                for(Int c = 0; c < uset_num; c++)
                {
                    gs_kuu[i][k][r][c] += ls_kuu[r][c];
                }
            }
            delete chol_sdd;
        }
    }
}
void GPBF::IOPITC_Sum(int t)
{
    MatDoub ls_kuu;
    VecDoub ls_zu;

    MatDoub D;
    D.resize(SLICE_SIZE, 3);

    for (int k = 0; k < field_num; k++)
    {
        for (int i = 0; i < par_num; i++)
        {
            for (int j = 0; j < SLICE_SIZE; j++)      // summarize data from t - SLICE_SIZE to t - 1
            {
                //init train data
                D[j][0] = ldata[i][k][j][0];
                D[j][1] = ldata[i][k][j][1];
                D[j][2] = obs[j+t-SLICE_SIZE][k];
            }

            //local summary in block i

            Cholesky *chol_sdd = gp[k].chol_pcov(D, SLICE_SIZE, uset, uset_num, chol_kuu[k]);
            gp[k].pitc_prep(D, SLICE_SIZE, chol_sdd, uset, uset_num, chol_kuu[k], ls_zu, ls_kuu);

            //global summary
            for(Int r = 0; r < uset_num; r++)
            {
                gs_zu[i][k][r] += ls_zu[r];
                for(Int c = 0; c < uset_num; c++)
                {
                    gs_kuu[i][k][r][c] += ls_kuu[r][c];
                }
            }
            delete chol_sdd;

            delete chol_suu[i][k];
            chol_suu[i][k] = gp[k].chol_cov(gs_kuu[i][k]);

            int local_num = 0;
            VecDoub tpose(3), ppose(3);
            for (int r = t; r < t + SLICE_SIZE && r <= step_num; r++)
            {
            	if (r == t)
					tpose = motionModel(bel[r-1][i], odom[r-1], odom[r]);
				else
					tpose = motionModel(ppose, odom[r-1], odom[r]);
            	ldata[i][k][local_num][0] = tpose[0];
	            ldata[i][k][local_num][1] = tpose[1];
	            ldata[i][k][local_num][2] = obs[r][k];
	            ppose = tpose;
            	local_num++;
            }
        	gp[k].pitc_regr_low3(uset, uset_num, gs_kuu[i][k], gs_zu[i][k], ldata[i][k], local_num, local_num, lmean[i][k], lvar);
        }
    }
    ldata_cnt = 1;
}
void GPBF::IOPITC_Regr(int t)
{
	for (int k = 0; k < field_num; k++)
        for (int i = 0; i < par_num; i++)
			gp[k].pitc_regr_low3(uset, uset_num, gs_kuu[i][k], gs_zu[i][k], test_data, par_num, par_num, pmean[i][k], pvar[i][k]);
    // all the data D has been summarized into global summary gs_zu and gs_kuu
    for (int k = 0; k < field_num; k++)
    {
        for (int i = 0; i < par_num; i++)
        {
            // prepare for the test data
            MatDoub tmp_cov, covx;

            //get predicive cov, pcov is \Sigma_XY, covx is \Sigma_XX
            gp[k].pitc_pcov_spec(uset, uset_num, chol_kuu[k], chol_suu[i][k], ldata[i][k], ldata_cnt, test_data, par_num, tmp_cov);
            gp[k].pitc_pcov(uset, uset_num, chol_kuu[k], chol_suu[i][k], ldata[i][k], ldata_cnt, covx);

            /* get the real pmean and pvar */
            Cholesky *chol_pxx = gp[k].chol_cov(covx);
            VecDoub zv(ldata_cnt), alpha(ldata_cnt), beta(ldata_cnt), covxy(ldata_cnt);

            for (int j = 0; j < ldata_cnt; j++)
            {
                zv[j] = ldata[i][k][j][2] - lmean[i][k][j];
            }

            /* alpha = \Sigma_XX * (z_X - \mu{X | D}) */
            chol_pxx->solve(zv, alpha);

            for (int j = 0; j < par_num; j++)
            {
                double tmpv = 0;
                for (int pt = 0; pt < ldata_cnt; pt++)
                	tmpv += alpha[pt] * tmp_cov[pt][j];

                pmean[i][k][j] += tmpv;

                /* beta = \Sigma_XX * \Sigma_{XY | D} */
                for (int pt = 0; pt < ldata_cnt; pt++)
                    covxy[pt] = tmp_cov[pt][j];
                chol_pxx->solve(covxy, beta);
                tmpv = 0;
                for (int pt = 0; pt < ldata_cnt; pt++)
                    tmpv += beta[pt] * tmp_cov[pt][j];

                pvar[i][k][j] -= tmpv;
            }

            delete chol_pxx;
        }
    }

    ldata_cnt++;
}
void GPBF::onlinePITC(VecDoub &weight, int t, int flag)
{
	for (int i = 0; i < par_num; i++)
    {
    	test_data[i][0] = bel_bar[i][0];
    	test_data[i][1] = bel_bar[i][1];
    }

	if (flag == 0)
	{
		for (int k = 0; k < field_num; k++)
		{
			for (int i = 0; i < par_num; i++)
        	{
        		// printf("samle #%d\n", i);
        		gp[k].pitc_regr_low3(uset, uset_num, gs_kuu[i][k], gs_zu[i][k], test_data, par_num, par_num, pmean[i][k], pvar[i][k]);
        	}
		}

		if (t % SLICE_SIZE == 0)    //time to update the global summary
			GlobalSum(t);
	}
	else if (t > 1)
	{
		if ((t - 1) % SLICE_SIZE == 0)    //time to update the global summary
		{
		    IOPITC_Sum(t);
		    for (int k = 0; k < field_num; k++)
                for (int i = 0; i < par_num; i++)
                    gp[k].pitc_regr_low3(uset, uset_num, gs_kuu[i][k], gs_zu[i][k], test_data, par_num, par_num, pmean[i][k], pvar[i][k]);
		}
		else
			IOPITC_Regr(t);
	}

	getObs(weight, t);
}

/*********************
 * batch PITC
 *********************/
void GPBF::batchPITC(VecDoub &weight, int t)
{
	initPITC();
	for (int i = 0; i < par_num; i++)
    {
    	test_data[i][0] = bel_bar[i][0];
    	test_data[i][1] = bel_bar[i][1];
    }

	for (int k = 0; k < field_num; k++)
    	for (int i = 0; i < par_num; i++)
			gp[k].pitc_regr_low3(uset, uset_num, gs_kuu[i][k], gs_zu[i][k], test_data, par_num, par_num, pmean[i][k], pvar[i][k]);
	for (int j = 1; j < t; j++)
		if (j % SLICE_SIZE == 0)    //time to update the global summary
			GlobalSum(j);
			
	for (int k = 0; k < field_num; k++)
		    	for (int i = 0; i < par_num; i++)
					gp[k].pitc_regr_low3(uset, uset_num, gs_kuu[i][k], gs_zu[i][k], test_data, par_num, par_num, pmean[i][k], pvar[i][k]);
	getObs(weight, t);
	if ((t - 1) % SLICE_SIZE == 0)
		sampleTraj(t);
	releaseMem();
}

static double diffclock(clock_t clock1,clock_t clock2)
{
    double diffticks=clock1-clock2;
    double diffms=(diffticks)/(CLOCKS_PER_SEC/1000);
    return diffms;
}

void GPBF::pFilter(int algo_id)
{
	double r, c, avg_x, avg_y, var_x, var_y;
	int offset;
	VecDoub weight;
	char algo_name[ALGO_NUM][STR_LEN] = {"trunc", "SoD", "OPITC", "IOPITC", "FGP", "PITC"};
	char str[STR_LEN];
	clock_t tstart, tend;

	FILE * fp_time;
	FILE * fp_err;
	FILE * fp_var;

	sprintf(str, "%s_%s_time", output_name, algo_name[algo_id]);
	fp_time = fopen(str, "w");
	sprintf(str, "%s_%s_err", output_name, algo_name[algo_id]);
	fp_err = fopen(str, "w");
	sprintf(str, "%s_%s_var", output_name, algo_name[algo_id]);
	fp_var = fopen(str, "w");

	weight.resize(par_num);
	for (int i = 1; i <= step_num; i++)
	{
		// FILE * fpb = fopen("bel_before", "w");
		// FILE * fpa = fopen("bel_after", "w");
		for (int m = 0; m < par_num; m++)
		{
			bel_bar[m] = motionModel(bel[i-1][m], odom[i-1], odom[i]);
			// fprintf(fpb, "%lf %lf\n", bel[i-1][m][0], bel[i-1][m][1]);
			// fprintf(fpa, "%lf %lf\n", bel_bar[m][0], bel_bar[m][1]);
			bel[i][m] = bel_bar[m];
		}

		tstart = clock();
		if (algo_id == 0)
			trunc(weight, i);
		else if (algo_id == 1)
			SoD(weight, i);
		else if (algo_id == 4)
			FGP(weight, i);
		else if (algo_id == 2)
			onlinePITC(weight, i, 0);
		else if (algo_id == 3)
			onlinePITC(weight, i, 1);
		else
			batchPITC(weight, i);
		tend = clock();

		// normalize the weight
		r = 0;
		for (int m = 0; m < par_num; m++)
            r += weight[m];
        for (int m = 0; m < par_num; m++)
            weight[m] /= r;

		//low_variance_sampler
		Ran uni(rand());
		r = uni.doub() / par_num;
		offset = 0;
		c = weight[0];
		for (int m = 0; m < par_num; m++)
		{
			double U = r + m * 1.0 / par_num;
			while (U > c)
			{
				offset++;
				c += weight[offset];
			}
			bel[i][m] = bel_bar[offset];
		}

		//output the result of current step
		avg_x = 0;
		avg_y = 0;
		var_x = 0;
		var_y = 0;
		for (int m = 0; m < par_num; m++)
		{
			avg_x += bel[i][m][0];
			avg_y += bel[i][m][1];
		}
		avg_x = avg_x / par_num;
		avg_y = avg_y / par_num;
		for (int m = 0; m < par_num; m++)
		{
			var_x += pow(avg_x - bel[i][m][0], 2);
			var_y += pow(avg_y - bel[i][m][1], 2);
		}
		var_x = sqrt(var_x / par_num);
		var_y = sqrt(var_y / par_num);
		fprintf(fp_time, "%d\n", (int)diffclock(tend, tstart));
		fprintf(fp_err, "%lf\n", sqrt(pow(real_loc[i][0] - avg_x, 2) + pow(real_loc[i][1] - avg_y, 2)));
		fprintf(fp_var, "%lf %lf\n", var_x, var_y);
		//fprintf(fp_var, "%lf %lf\n", avg_x, avg_y);
		//debug information
		printf("Step #%d:\n", i);
		// printf("Current location: (%lf, %lf)\n", real_loc[i][0], real_loc[i][1]);
		// printf("Estimated location: (%lf, %lf)\n", avg_x, avg_y);
		// printf("Estimation variance: %lf, %lf\n", var_x, var_y);
		// printf("Distance: %lf\n*****\n", sqrt(pow(real_loc[i][0] - avg_x, 2) + pow(real_loc[i][1] - avg_y, 2)));
	}
	fclose(fp_err);
	fclose(fp_var);
	fclose(fp_time);
}

void GPBF::deadReckon()
{
	FILE * fp_err;
	FILE * fp_var;
	char str[STR_LEN];
	double avg_x, avg_y, var_x, var_y;

	sprintf(str, "%s_dr_err", output_name);
	fp_err = fopen(str, "w");
	sprintf(str, "%s_dr_var", output_name);
	fp_var = fopen(str, "w");
	sprintf(str, "%s_dm_err", output_name);
	FILE * fp_derr = fopen(str, "w");

	for (int i = 1; i <= step_num; i++)
	{
		for (int m = 0; m < par_num; m++)
		{
			bel_bar[m] = motionModel(bel[i-1][m], odom[i-1], odom[i]);
			bel[i][m] = bel_bar[m];
		}

		//output the result of current step
		avg_x = 0;
		avg_y = 0;
		var_x = 0;
		var_y = 0;
		for (int m = 0; m < par_num; m++)
		{
			avg_x += bel[i][m][0];
			avg_y += bel[i][m][1];
		}
		avg_x = avg_x / par_num;
		avg_y = avg_y / par_num;
		for (int m = 0; m < par_num; m++)
		{
			var_x += pow(avg_x - bel[i][m][0], 2);
			var_y += pow(avg_y - bel[i][m][1], 2);
		}
		var_x = sqrt(var_x / par_num);
		var_y = sqrt(var_y / par_num);
		fprintf(fp_err, "%lf\n", sqrt(pow(real_loc[i][0] - avg_x, 2) + pow(real_loc[i][1] - avg_y, 2)));
		fprintf(fp_derr, "%lf\n", sqrt(pow(real_loc[i][0] - es_loc[i][0], 2) + pow(real_loc[i][1] - es_loc[i][1], 2)));
		fprintf(fp_var, "%lf %lf\n", var_x, var_y);
	}

	fclose(fp_err);
	fclose(fp_var);
}

void GPBF::releaseMem()
{
	for (int k = 0; k < field_num; k++)
	{
		delete chol_kuu[k];
		for (int i = 0; i < par_num; i++)
			delete chol_suu[i][k];
	}
}

void GPBF::filter()
{
	//initial pose
	for (int i = 0; i < par_num; i++)
	    bel[0][i] = randomPose(real_loc[0]);

	deadReckon();

	for (int i = 0; i < par_num; i++)
		for (int j = 0; j < field_num; j++)
		{
			pmean[i][j].resize(par_num);
			pvar[i][j].resize(par_num);
		}

	//start filtering
	for (int i = 0; i < algo_num; i++)
	{
		if (algo[i] == 0)
		{
			pFilter(0);
			printf("truncation completed.\n----------------------------\n");
		}
		else if (algo[i] == 4)
		{
			pFilter(4);
			printf("Full GP completed.\n----------------------------\n");
		}
		else if (algo[i] == 1)
		{
			initSoD();
			pFilter(1);
			printf("SoD completed.\n----------------------------\n");
		}
		else if (algo[i] == 2)
		{
			initPITC();
			pFilter(2);
			releaseMem();
			printf("online PITC completed.\n----------------------------\n");
		}
		else if (algo[i] == 3)
		{
			initPITC();
			pFilter(3);
			releaseMem();
			printf("improved online PITC completed.\n----------------------------\n");
		}
		else
		{
			pFilter(5);
			printf("batch PITC completed.\n----------------------------\n");
		}
	}
}
