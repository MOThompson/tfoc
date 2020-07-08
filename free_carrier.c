/* Routines for calculating the k constant and alpha factor for
 * free-electron absorption in Si */

/* ------------------------------ */
/* Standard include files         */
/* ------------------------------ */
#define _POSIX_SOURCE						/* Always require POSIX standard */
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>

/* ------------------------------ */
/* Local include files            */
/* ------------------------------ */
#define TFOC_CODE
#include "tfoc.h"
#include "gcc_help.h"

/* ------------------------------- */
/* My local typedef's and defines  */
/* ------------------------------- */

/* ------------------------------- */
/* My external function prototypes */
/* ------------------------------- */

/* ------------------------------- */
/* My internal function prototypes */
/* ------------------------------- */

/* ------------------------------- */
/* My usage of other external fncs */
/* ------------------------------- */

/* ------------------------------- */
/* My share of global externals    */
/* ------------------------------- */

/* ------------------------------- */
/* Locally defined global vars     */
/* ------------------------------- */

/** ======================================================
-- Blind application of the free carrier rules end up
-- with absorption coefficients that imply a completely
-- reflecting Si surface as the temperature increases to
-- near melting.  Given that this isn't seen experimentally,
-- the absorption coefficient is artifically "capped" at a
-- value of 3E4 (300 nm) giving a reflectance near Brewster's
-- angle of 10-11%.  Needs serious experimental measurements!
====================================================== **/
#define	MAX_ALPHA	(3E4)

/* ===========================================================================
-- Mobility model from
--   Electron and Hole Mobility in Silicon at Large Operating
--   Temperatures - Part I: Bulk Mobility, Susanna Reggiani, Marina
--   Valdinoci, Luigi Colalongo, Massimo Rudan, Giorgio Baccarani, 
--   Andreas D. Stricker, Fridolin Illien, Norbert Felber, Wolfgang 
--   Fichtner and Lucia Zullino, IEEE Transactions on Electron Devices
--   Vol 49, No. 3, March 2002.
=========================================================================== */
typedef struct _MOBILITY_PARMS {
	char *name;
	double u_max;
	double c;
	double gamma;
	double u0d_c, u0d_p;				/* u0d = u0d_c * (T/300)^u0d_p	*/
	double u0a_c, u0a_p;
	double u1d_c, u1d_p;
	double u1a_c, u1a_p;
	double Cr1_c, Cr1_p;				/* Cr1 = Cr1_c * (T/300)^Cr1_p	*/
	double Cr2_c, Cr2_p;
	double Cs1_c, Cs1_p;
	double Cs2_c, Cs2_p;
	double alpha_1, alpha_2;
} MOBILITY_PARMS;

MOBILITY_PARMS As = {"Arsenic",
							1441, 0.07, 2.45,			/* u_max, c, gamma, */
							55.0,		-0.6,				/* u0d */
							132.0,	-1.3,				/* u0a */
							42.4,		-0.5,				/* u1d */
							73.5,		-1.25,			/* u1a */
							8.90E16,	 3.65,			/* Cr1 */
							1.22E17,	 2.65,			/* Cr2 */
							2.9E20,	 0.0,				/* Cs1 */
							7.0E20,	 0.0,				/* Cs2 */
							0.68, 0.72 };				/* alpha_1, alpha_2 */
MOBILITY_PARMS P  = {"Phosphorus",
							1441, 0.07, 2.45,			/* u_max, c, gamma, */
							62.2,		-0.7,				/* u0d */
							132.0,	-1.3,				/* u0a */
							48.6,		-0.7,				/* u1d */
							73.5,		-1.25,			/* u1a */
							8.50E16,	 3.65,			/* Cr1 */
							1.22E17,	 2.65,			/* Cr2 */
							4.0E20,	 0.0,				/* Cs1 */
							7.0E20,	 0.0,				/* Cs2 */
							0.68, 0.72 };				/* alpha_1, alpha_2 */
MOBILITY_PARMS B  = {"Boron",
							470.5, 0.0, 2.16,			/* u_max, c, gamma, */
							90.0,		-1.3,				/* u0d */
							44.0,		-0.7,				/* u0a */
							28.2,		-2.0,				/* u1d */
							28.2,		-0.8,				/* u1a */
							1.30E18,	 2.2,				/* Cr1 */
							2.45E17,	 3.1,				/* Cr2 */
							1.1E18,	 6.2,				/* Cs1 */
							6.1E20,	 0.0,				/* Cs2 */
							0.77, 0.719 };				/* alpha_1, alpha_2 */

#define	MU(M,ND,NA,T)	(u0(M,ND,NA,T)+(uL(M,T)-u0(M,ND,NA,T))/(1+pow(ND/Cr1(M,T),M->alpha_1)+pow(NA/Cr2(M,T),M->alpha_2)) - u1(M,ND,NA,T)/(1+pow(ND/Cs1(M,T)+NA/Cs2(M,T),-2)))
#define	uL(M,T)			M->u_max*pow(T/300,-M->gamma+M->c*T/300)
#define	u0d(M,T)			(M->u0d_c*pow(T/300,M->u0d_p))
#define	u0a(M,T)			(M->u0a_c*pow(T/300,M->u0a_p))
#define	u1d(M,T)			(M->u1d_c*pow(T/300,M->u1d_p))
#define	u1a(M,T)			(M->u1a_c*pow(T/300,M->u1a_p))
#define	Cr1(M,T)			( (M->Cr1_p==0) ? M->Cr1_c : (M->Cr1_c*pow(T/300,M->Cr1_p)) )
#define	Cr2(M,T)			( (M->Cr2_p==0) ? M->Cr2_c : (M->Cr2_c*pow(T/300,M->Cr2_p)) )
#define	Cs1(M,T)			( (M->Cs1_p==0) ? M->Cs1_c : (M->Cs1_c*pow(T/300,M->Cs1_p)) )
#define	Cs2(M,T)			( (M->Cs2_p==0) ? M->Cs2_c : (M->Cs2_c*pow(T/300,M->Cs2_p)) )

#define	u0(M,ND,NA,T)	(u0d(M,T)*ND+u0a(M,T)*NA)/(ND+NA+1)
#define	u1(M,ND,NA,T)	(u1d(M,T)*ND+u1a(M,T)*NA)/(ND+NA+1)

/* ===========================================================================
-- Free carrier absorption from Dieter K. Schroder, Semiconductor Material and
-- Device Characterization, Wiley and Sons, 1990, p. 83.
=========================================================================== */
#define	kBoltz	8.62E-5					/* Boltzmann's constant eV/K			*/
#define	NC300		2.8E19					/* Effective density at 300K			*/
#define	NV300		1.04E19					/* Effective density at 300K			*/
#define	EG300		1.08						/* Bandgap at 300K						*/
#define	EG0		1.1255					/* Bandgap at 0K (calc from alpha/beta) */
#define	EGALPHA	4.73E-04					/* Constant in bandgap formula		*/
#define	EGBETA	636						/* Constant in bandgap formula		*/
#define	Eg(T)		(EG0-EGALPHA*pow(T,2)/(EGBETA+(T)))
#define	NI2(T)	(NC300*NV300)*pow((T)/300,3)*exp(-Eg(T)/kBoltz/(T))

#define	n_Si		3.412						/* Index of Si over range 4-12 um	*/

/* ===========================================================================
-- Effective mass and degeneracy of electrons and holes.
-- Electrons is relatively well established, but hole data is basically
-- a good guess currently.  Matches "effective" number used in all
-- bloody text books.
=========================================================================== */
#define	N_MSTAR_L		(0.19)			/* Longitudinal effective mass		*/
#define	N_MSTAR_T		(0.98)			/* Transverse effective mass			*/
#define	P_MSTAR_H		(0.50)			/* Heavy hole effective mass			*/
#define	P_MSTAR_L		(0.16)			/* Light hole effective mass			*/
#define	N_DEGENERACY_L	(4)				/* Degeneracy of logitudinal e-		*/
#define	N_DEGENERACY_T	(2)				/* Degeneracy of transverse e-		*/
#define	P_DEGENERACY_H	(4)				/* Degeneracy of heavy holes			*/
#define	P_DEGENERACY_L	(2)				/* Degeneracy of light holes			*/

#define	N_MSTAR_AVG_2		(sqrt(1/((N_DEGENERACY_L/pow(N_MSTAR_L,2)+N_DEGENERACY_T/pow(N_MSTAR_T,2))/(N_DEGENERACY_L+N_DEGENERACY_T))))
#define	P_MSTAR_AVG_2		(sqrt(1/((P_DEGENERACY_L/pow(P_MSTAR_L,2)+P_DEGENERACY_H/pow(P_MSTAR_H,2))/(P_DEGENERACY_L+P_DEGENERACY_H))))
#define	N_MSTAR_AVG_1		(1/((N_DEGENERACY_L/N_MSTAR_L+N_DEGENERACY_T/N_MSTAR_T)/(N_DEGENERACY_L+N_DEGENERACY_T)))
#define	P_MSTAR_AVG_1		(1/((P_DEGENERACY_L/P_MSTAR_L+P_DEGENERACY_H/P_MSTAR_H)/(P_DEGENERACY_L+P_DEGENERACY_H)))

#define	N_MSTAR_EFF	(0.26)				/* Electron effective mass				*/
#define	P_MSTAR_EFF	(0.50)				/* Hole effective mass					*/

/**
 * Modified 8/9/2006
 *   Previous version used p_mstar_avg_2 for the hole mass as well as
 *   electron.  After discussion, conclude there is a significant
 *   difference between the two bands in the holes and just the two
 *   axes of the single electron band.  May make a significant
 *   difference.
    setv n_mstar_l = 0.19
    setv n_mstar_t = 0.98
    setv p_mstar_h = 0.50
    setv p_mstar_l = 0.16
    setv n_degeneracy_t = 2
    setv n_degeneracy_l = 4
    setv p_degeneracy_h = 4
    setv p_degeneracy_l = 2
    eval (sqrt(1/((n_degeneracy_l/pow(n_mstar_l,2)+n_degeneracy_t/pow(n_mstar_t,2))/(n_degeneracy_l+n_degeneracy_t))))
       0.23054515
    eval (sqrt(1/((p_degeneracy_l/pow(p_mstar_l,2)+p_degeneracy_h/pow(p_mstar_h,2))/(p_degeneracy_l+p_degeneracy_h))))
       0.25247776
    eval (1/((n_degeneracy_l/n_mstar_l+n_degeneracy_t/n_mstar_t)/(n_degeneracy_l+n_degeneracy_t)))
       0.25981395
    eval (1/((p_degeneracy_l/p_mstar_l+p_degeneracy_h/p_mstar_h)/(p_degeneracy_l+p_degeneracy_h)))
       0.29268292
 *   **/

/* Choose a model -- routine below can modify */
#if 1							/* Model using both bands for holes - seems closer to expt data */
  #define	P_MSTAR	(P_MSTAR_AVG_2)		/* For holes, let both bands contribute equally */
  #define	N_MSTAR	(N_MSTAR_AVG_2)		/* For electrons, both long and traverse modes of a single band contribute */
#else							/* Model using heavy holes only - as the band most occupied by holes */
  #define	P_MSTAR	(P_MSTAR_EFF)			/* For holes, correct is just the heavy *band* */
  #define	N_MSTAR	(N_MSTAR_AVG_2)		/* For electrons, both long and traverse modes of a single band contribute */
#endif

#define	PI			(3.14159265359)		/* Constants!!!							*/

/* Constants from carrier-carrier scattering mobility model */
#define mun_min	(55.24)
#define mun_max	(1429.23)
#define mup_min	(49.70)
#define mup_max	(379.37)
#define nrefn		(1.072E17)
#define nrefp		(1.606E17)
#define nun			(-2.3)
#define xin			(-3.8)
#define alphan		(0.73)
#define nup			(-2.2)
#define xip			(-3.7)
#define alphap		(0.70)

/* ===========================================================================
-- Concentration dependent mobility data
-- Existing data from table in Avante TSUPREM manuals
=========================================================================== */
static struct {
	double conc, mu_n, mu_p;
} raw[] = {
	{1.0E14,	1350,	495},		{2.0E14,	1345,	495},		{4.0E14,	1335,	495},		{6.0E14,	1320,	495},		{8.0E14,	1310,	495},
	{1.0E15,	1300,	491.1},	{2.0E15,	1248,	487.3},	{4.0E15,	1200,	480.1},	{6.0E15,	1156,	473.3},	{8.0E15,	1115,	466.9},
	{1.0E16,	1076,	460.9},	{2.0E16,	960,	434.8},	{4.0E16,	845,	396.5},	{6.0E16,	760,	369.2},	{8.0E16,	720,	348.3},
	{1.0E17,	675,	331.5},	{2.0E17,	524,	279.0},	{4.0E17,	385,	229.8},	{6.0E17,	321,	203.8},	{8.0E17,	279,	186.9},
	{1.0E18,	252,	178},		{2.0E18,	182.5,130.0},	{4.0E18,	140.6,90.0},	{6.0E18,	113.6,74.5},	{8.0E18,	99.5,	66.6},
	{1.0E19,	90.5,	61.0},	{2.0E19,	86.9,	55.0},	{4.0E19,	83.4,	53.7},	{6.0E19,	78.8,	52.9},	{8.0E19,	71.6,	52.4},
	{1.0E20,	67.8,	52.0},	{2.0E20,	52.0,	50.8},	{4.0E20,	35.5,	49.6},	{6.0E20,	23.6,	48.9},	{8.0E20,	19.0,	48.4},
	{1.0E21,	17.8,	48.0}
};
#define	MU_NPT	(sizeof(raw)/sizeof(*raw))

#define	USE_FULL_FE_MODEL								/* Full model for free-electron */
/* #define USE_SIMPLE_CODE */

/* ===========================================================================
-- Routine to modify the model for the effective mass in the free
-- carrier absorption calculations.
--
-- Usage: void fc_set_mstar_mode(int mode);
--
-- Inputs: mode - model to use for the effective mass
--           0 ==> Use default defined in the code
--           1 ==> Use all bands of electrons and holes (default)
--           2 ==> Use all bands for electrons but only heavy holes
=========================================================================== */
static double p_mstar = -1;			/* For use below - local constant */
static double n_mstar = -1;

void fc_set_mstar_mode(int mode) {

	switch (mode) {
		case 1:
			p_mstar = P_MSTAR_AVG_2;		/* For holes, use both holes in averaged mode */
			n_mstar = N_MSTAR_AVG_2;		/* For electrons, use both holes in averaged mode */
			break;
		case 2:
			p_mstar = P_MSTAR_EFF;			/* For holes, use just heavy band since it will occupied */
			n_mstar = N_MSTAR_AVG_2;		/* For electrons, both long and traverse modes of a single band contribute */
			break;
		default:
			p_mstar = P_MSTAR;				/* Anything else, use program default */
			n_mstar = N_MSTAR;
	}
	return;
}

/* ===========================================================================
-- Routines to return the mobility as a function of dopant concentrations
-- and temperature.  Can use any of three models:
--     KLAASSEN - full parameterization over temperature and doping
--     SPLINE   - spline fit to room temperature data from Silvaco
--     SIMPLE   - simplified model as function of doping and temperature
--
-- Usage: double mu_n(double Nd, double Na, double T, FC_MODE mode);
--        double mu_p(double Nd, double Na, double T, FC_MODE mode);
--
-- Input: Nd - donor concentration (ionized impurities)			/cm^3
--        Na - acceptor concentration (ionized impurities)		/cm^3
--        T  - absolute temperature										K
--        mode - type of mobility data (see fc_alpha, fc_k decl)
--
-- Output: none
--
-- Return: Mobility estimate in cm^2/V-s
=========================================================================== */
double mu_n(double nd, double na, double T, FC_MODE mode) {
	double mu;
	static void *mu_n_spline=NULL;

	if (nd <= 0) nd = 1;									/* Avoid blowup */
	if (na <= 0) na = 1;
	switch (mode) {
		case KLAASSEN_MU:
			mu = MU((&As),nd,na,T);
			break;
		case SPLINE_MU:
			if (mu_n_spline == NULL) {
				double x[MU_NPT], y[MU_NPT];
				int i;
				for (i=0; i<MU_NPT; i++) x[i] = log(raw[i].conc)/log(10.0);
				for (i=0; i<MU_NPT; i++) y[i] = raw[i].mu_n;
				mu_n_spline = GVFitSpline(NULL, x,y, MU_NPT, 0);
			}
			mu = GVEvalSpline(mu_n_spline, log(nd)/log(10.0));
			break;
		default:
			mu = mun_min+(mun_max*pow(T/300.0,nun)-mun_min)/(1+pow(T/300,xin)*pow(nd/nrefn,alphan));
	}
	return mu;
}

double mu_p(double nd, double na, double T, FC_MODE mode) {
	double mu;
	static void *mu_p_spline=NULL;

	if (nd <= 0) nd = 1;									/* Avoid blowup */
	if (na <= 0) na = 1;
	switch (mode) {
		case KLAASSEN_MU:
			mu = MU((&B),nd,na,T);
			break;
		case SPLINE_MU:
			if (mu_p_spline == NULL) {
				double x[MU_NPT], y[MU_NPT];
				int i;
				for (i=0; i<MU_NPT; i++) x[i] = log(raw[i].conc)/log(10.0);
				for (i=0; i<MU_NPT; i++) y[i] = raw[i].mu_p;
				mu_p_spline = GVFitSpline(NULL, x,y, MU_NPT, 0);
			}
			mu = GVEvalSpline(mu_p_spline, log(na)/log(10.0));
			break;
		default:
			mu = mup_min+(mup_max*pow(T/300.0,nup)-mup_min)/(1+pow(T/300,xip)*pow(na/nrefp,alphap));
	}
	return mu;
}


/* ===========================================================================
-- Usage: fc_alpha(double doping, double excess, double T, double lambda, FC_MODE mode);
--        fc_k(double doping, double excess, double T, double lambda, FC_MODE mode);
--
-- Inputs: doping - doping concentration.  - => n-type, + => p-type
--                    +1E19 ==> p type 10^19 concentration
--                    -1E15 ==> n type 10^14 concentration
--                  Intrinsic is added to the extrinsic (fnc of T)
--                  All dopants are considered as fully ionized
--         excess - excess n/p carriers (as by optical generation for example)
--                  Added to n/p calculated from doping, equiv to more intrinsic
--         T      - absolute temperature (K)
--         lambda - wavelength in microns
--         mode   - One of the defined constants FIXED_K     = 0
--                                               SIMPLE_MU   = 1
--                                               SPLINE_MU   = 2
--                                               KLAASSEN_MU = 3 <== RECOMMENDED
--
-- Output: none
--
-- Return: Free-carrier absoprtion length (cm-1) or imaginary component (k).
--            alpha^-1 = lambda/(4*pi*k)
--
-- Modification: Limit the carrier concentration to 10% - 5E21
=========================================================================== */
double fc_alpha(double doping, double excess, double T, double lambda, FC_MODE mode) {

	double N_D, N_A;
	double n_n, n_p, alpha, alpha_n, alpha_p;

/* On first call, set the MSTAR mode if not already set */
	if (p_mstar <= 0.0) fc_set_mstar_mode(0);

	if (T > 1683) T = 1683;											/* Limit the temperature */
	if (T < 77)   T = 77;

/* Calculate the total number of free carriers */
	if (doping < 0) {													/* n-type						*/
		N_D = fabs(doping);	
		N_A = 1;
		n_n = N_D + sqrt(NI2(T));
		n_p = NI2(T)/n_n;
	} else {																/* p-type						*/
		N_A = fabs(doping);	
		N_D = 1;
		n_p = N_A + sqrt(NI2(T));
		n_n = NI2(T)/n_p;
	}
	n_n += excess; n_p += excess;									/* Excess carrier pairs present */
	if (n_p > 5E21) n_p = 5E21;									/* At these levels, ni^2 doesn't matter */
	if (n_n > 5E21) n_n = 5E21;

#ifdef USE_FULL_FE_MODEL
	alpha_p = 5.27E-17 / (n_Si*mu_p(N_D,N_A,T,mode)*pow(p_mstar,2)) * n_p * pow(lambda,2);
	alpha_n = 5.27E-17 / (n_Si*mu_n(N_D,N_A,T,mode)*pow(n_mstar,2)) * n_n * pow(lambda,2);
	alpha = alpha_n + alpha_p;

#elif defined USE_SIMPLE_CODE
	alpha = 1E-18*(n_n+2.7*n_p)*pow(lambda,2);				/* In cm^{-1} */
#endif

	if (alpha > MAX_ALPHA) alpha = MAX_ALPHA;					/* Limit absorption length to 100 nm */
	return alpha;
}

double fc_k(double doping, double excess, double T, double lambda, FC_MODE mode) {
	return fc_alpha(doping, excess, T,lambda,mode) * (lambda*1E-4) / (4*PI);		/* Convert to k value (from cm-1) */
}

/* =========================================================================== */
/* =========================================================================== */
/* =========================================================================== */
/* =========================================================================== */
/* =========================================================================== */

#ifdef TEST

int main(int argc, char *argv[]) {

	double lnd,nd,T;

	FILE *funit;

	funit = fopen("fig4.dat", "w");
	fprintf(funit, "/* Log(N)\tAs\tP\tB\n");
	for (lnd=14; lnd<=22; lnd+=0.01) {
		nd = pow(10.0,lnd);
		fprintf(funit, "%g\t%g\t%g\t%g\n", lnd, MU((&As),nd,1,300), MU((&P),nd,1,300), MU((&B),1,nd,300));
	}
	fclose(funit);

	funit = fopen("fig5_As.dat", "w");
	fprintf(funit, "/* T\t9.0E15\t2.0E17\t1.0E18\t2.5E18\n");
	for (T=200; T<=1685; T+=1) {
		fprintf(funit, "%g\t%g\t%g\t%g\t%g\n", T, MU((&As),9.0E15,1,T), MU((&As),2.0E17,1,T), MU((&As),1.0E18,1,T), MU((&As),2.5E18,1,T));
	}
	fclose(funit);

	funit = fopen("fig5_P.dat", "w");
	fprintf(funit, "/* T\t9.0E15\t2.0E17\t1.0E18\t2.5E18\n");
	for (T=200; T<=1685; T+=1) {
		fprintf(funit, "%g\t%g\t%g\t%g\t%g\n", T, MU((&P),9.0E15,1,T), MU((&P),2.0E17,1,T), MU((&P),1.0E18,1,T), MU((&P),2.5E18,1,T));
	}
	fclose(funit);

	funit = fopen("fig6.dat", "w");
	fprintf(funit, "/* T\5.60E16\t3.0E17\t7.0E17\t3.2E18\n");
	for (T=200; T<=1685; T+=1) {
		fprintf(funit, "%g\t%g\t%g\t%g\t%g\n", T, MU((&B),1,5.6E16,T), MU((&B),1,3.0E17,T), MU((&B),1,7.0E17,T), MU((&B),1,3.2E18,T));
	}
	fclose(funit);

	funit = fopen("fig7.dat", "w");
	fprintf(funit, "/* Log(N)\t300\t400\t500\t600\t700\t800\n");
	for (lnd=14; lnd<=22; lnd+=0.01) {
		nd = pow(10.0,lnd);
		fprintf(funit, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", lnd, MU((&As),1,nd,300.0), MU((&As),1,nd,400.0), MU((&As),1,nd,500.0), MU((&As),1,nd,600.0), MU((&As),1,nd,700.0), MU((&As),1,nd,800.0));
	}
	fclose(funit);

	funit = fopen("fig8.dat", "w");
	fprintf(funit, "/* Log(N)\t300\t400\t500\t600\t700\t800\n");
	for (lnd=14; lnd<=22; lnd+=0.01) {
		nd = pow(10.0,lnd);
		fprintf(funit, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", lnd, MU((&B),nd,1,300.0), MU((&B),nd,1,400.0), MU((&B),nd,1,500.0), MU((&B),nd,1,600.0), MU((&B),nd,1,700.0), MU((&B),nd,1,800.0));
	}
	fclose(funit);

	funit = fopen("fig9.dat", "w");
	fprintf(funit, "/* T\tn_i\tk\talpha (cm^-1)\tabsorption length (um)\n");
	for (T=200; T<=1685; T+=10) {
		fprintf(funit, "%g\t%g\t%g\t%g\t%g\n", T, sqrt(NI2(T)), fc_k(1.0,0,T,10.6,KLAASSEN_MU), fc_alpha(1.0,0,T,10.6,KLAASSEN_MU), 1E4/fc_alpha(1.0,T,10.6,KLAASSEN_MU));
	}
	fclose(funit);

	return(0);
}

int main1(int argc, char *argv[]) {

	MOBILITY_PARMS *model;
	double nd,na,T;

	if (argc < 5) {
		fprintf(stderr, "Usage: mu <ND> <NA> <T> [As | P | B]\n");
		return(0);
	}

	nd = atof(argv[1]);
	na = atof(argv[2]);
	T  = atof(argv[3]);
	if (stricmp(argv[4], "P") == 0)			{ model = &P; }
	else if (stricmp(argv[4], "As") == 0)	{ model = &As; }
	else if (stricmp(argv[4], "B") == 0)	{ model = &B; }
	else {
		fprintf(stderr, "You must specify either As, P or B as the 4th parameter\n");
		return(0);
	}

	printf("T=%g\tNd=%g\tNa=%g\tDoping=%s\tMobility = %g\n", T,nd,na,model->name,MU(model,nd,na,T));
	return(0);
}

#endif
