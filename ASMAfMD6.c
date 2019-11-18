/************************************************************/
/*                                                          */
/*      ASMAfMD6.c                                          */
/*                                                          */
/*      2019/10/3                                           */
/************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <signal.h>

#define N_assm    8       // cell assemblies 
#define N_T       19      // number of cells within cell assemlbies 
#define DT        1e-4    // discrete time 
#define delay     500     // delay between N_S and N_M
#define cmPY      500e-12 // membrame capacitance of Pyramidal cell 
#define cmSB      243e-12 // non-use 
#define cmLB      115e-12 // large basket cell 
#define cmG       10e-12  // non-use 
#define gmPY      25.0e-9 // membrame conductance of Pyramidal cell
#define gmSB      9.7e-9  // non-use 
#define gmLB      8.2e-9  //  
#define gmG       20.0e-9 // non-use  
#define gGap      20.0e-9 // non-use
#define gAMPA     0.5e-9  // maximal conductance for AMPA-receptors 
#define gGABA     0.7e-9  // for GABA-receptors 

#define gGABAb    1e-9    // non-use 
#define K1        9e4     // non-use  
#define K2        1.2     // non-use 
#define K3        180     // non-use 
#define K4        34      // non-use 
#define Kd        100     // non-use 
#define nBS       4       // non-use

#define UPYact    -0.01   // action potential of pyramidal cell
#define UPYres    -0.065  // reseting potential  
#define USBact    -0.01   // non-use
#define USBres    -0.07   //    
#define ULBact    -0.01   // of large basket cell
#define ULBres    -0.07   //  
#define UGres     -0.07   // non-use  

#define u_AMPA     0.0    // reversal potential for AMPA receptor
#define u_GABA    -0.08   // reversal potential for GABA_A receptor
#define u_GABAb   -0.095  // non-use 
#define steep_PY   240.0  // steepness of sigmoid function for N_S
#define thres_PY   -0.033 // threshold 
#define steep_PY2  240.0  // for N_M 
#define thres_PY2  -0.033 //  
#define steep_SB   260.0  // non-use 
#define thres_SB   -0.033 // non-use     
#define steep_SB2  300.0  // non-use 
#define thres_SB2  -0.035 // non-use     
#define steep_LB   300.0  // for N_S  
#define thres_LB   -0.031 //   
#define steep_LB2   300.0 // for N_M  
#define thres_LB2  -0.031 //   

#define alph_AMPA  1.1e6  // AMPA-channels 
#define beta_AMPA  190.0  // 
#define alph_GABA  5.0e6  // GABA-channels 
#define beta_GABA  180.0  // 
#define Glut_c     0.001  // intrasyn. glutamate conc. in N_S   
#define GABA_cS    0.001  // non-use     
#define GABA_cL1   0.001  // intrasyn. GABA conc. in N_S 
#define GABA_cL2   0.001  // in N_M

double m_G=0e9;           // non-use  
#define uG_trn     -0.07  // non-use

#define theta_inp0   0     
#define theta_inp1   1    
#define theta_inp2   2     
#define theta_inp3   3     
#define theta_inp4   4      
#define theta_inp5   5     
#define theta_inp6   6     
#define theta_inp7   7    
#define theta_inp8   8    

#define int_inp0_0   0e-12        
#define int_inp0_1   int_inp0_0  
#define int_inp0_2   int_inp0_0   
#define int_inp0_3   600e-12     
#define int_inp0_4   int_inp0_0   
#define int_inp0_5   int_inp0_0   
#define int_inp0_6   int_inp0_0   
#define int_inp0_7   int_inp0_0   
#define inp_prob     1.0        // input probability [0,1]

double delta_T=2e2;             // extrasynaptic GABA_A receptors in N_S
double delta_T2=2e2;            // in N_M
#define delta_SB     0.0        // non-use
#define delta_LB     0.0        // non-use

#define tau_area     1.0          
#define COLUMN       3           
#define NEURON       1           

int t;                            
#define OUT       onset_0-10000   
#define PAT       "ON"          
#define INTV      99              
#define PERIOD    OUT+30000       

#define onset_0   20000           
#define period_0  10000           

#define NEURONi      2      
#define NEURONj      4     
#define NEURONw      3     

#define wLPPV1		1.0		
#define wLPPV2		1.0	    

// amb. GABA conc.   
void dfsGABA_ext(void);            // diffusion equations  
double I_GABA[N_assm+2][N_T+2];    // curretn 
double GABA_ext[N_assm+2][N_T+2];  // local amb. GABA conc. 
double GABA_extw[N_assm+2][N_T+2]; //  
double gamma=10.0;                 // decay 

double GABA_c0=20E-07;             // basal amb. GABA conc. in N_S
double GABAamb_max=4E-06;          // non-use
double GABAamb_min=0E-07;	       // non-use
double GABA_V2=20E-07;             // in N_M

double rEXT[N_assm+2][N_T+2];      // fraction of opened extrasyn. GABA receptors in N_S
double difrEXT(int,int);         
double rEXT2[N_assm+2][N_T+2];     // in N_S 
double difrEXT2(int,int);        

// N_S
double I_MGB1[N_assm+2];           // input current 
double uPY1[N_assm+2][N_T+2];      // membrane pot.   
double uSB1[N_assm+2][N_T+2];      // non-use 
double uLB1[N_assm+2][N_T+2];   
double uG[N_assm+2][N_T+2];        // non-use
double vPY1[N_assm+2][N_T+2][PERIOD+2]; // action potential   
double vSB1[N_assm+2][N_T+2][PERIOD+2]; // non-use 
double vLB1[N_assm+2][N_T+2][PERIOD+2]; //  
double rPY1[N_assm+2][N_T+2];           // fraction of opened intrasyn. AMPA receptors
double rPY1d[N_assm+2][N_T+2][PERIOD+2];// delay to Nmot 
double rSB1[N_assm+2][N_T+2];           // non-use 
double sF[N_assm+2][N_T+2];             // non-use  
double rLB1[N_assm+2][N_T+2];           // for GABA-receptors
double GlutPY_c1[N_assm+2][N_T+2];      // synaptic glutamate release 
double GlutSB_c1[N_assm+2][N_T+2];      // non-use 
double GlutLB_c1[N_assm+2][N_T+2];      // synaptic GABA  
double duPY_leak1[N_assm+2][N_T+2];     // membrane pot. increase by leak current
double duPY_rec_1[N_assm+2][N_T+2];     // by recurrent excitatory current
double duPY_fed_1[N_assm+2][N_T+2];     // non-use
double duPY_lat_1[N_assm+2][N_T+2];     // by lateral inhibitory current
double duPY_topdown[N_assm+2][N_T+2];   // by topdown (N_M to N_S) excitatory current 
double duPY_ext_1[N_assm+2][N_T+2];     // by GABAergc tonic inhibitory current
double duPY_MGB1[N_assm+2];             // by external excitatory current
double I_leak1;                         // leak current
double I_rec_1;                         // recurrent excitatory current
double I_fed_1;                         // non-use
double I_lat_1;                         // lateral inhibitory current
double I_topdown;                       // topdown (N_M to N_S) excitatory current 
double I_ext_1;                         // GABAergc tonic inhibitory current
double I_MG1;                           // external excitatory current

double duSB_leak1[N_assm+2][N_T+2];     // non-use
double duSB_1[N_assm+2][N_T+2];         // non-use
double duSB_topdown[N_assm+2][N_T+2];   // non-use 
double duSB_ext_1;                      // non-use
double I_leakSB1;                       // non-use
double I_fed_SB1;                       // non-use
double I_ext_SB1;                       // non-use

double duLB_leak1[N_assm+2][N_T+2];     
double duLB1[N_assm+2][N_T+2];
double duLB1_topdown[N_assm+2][N_T+2];   
double duLB_ext1;                       
double I_leakLB1;                       
double I_fed_LB1;                       
double I_ext_LB1;                       

double duG_leak[N_assm+2][N_T+2];       // non-use
double duG_P[N_assm+2][N_T+2];          // non-use
double duG_Ia[N_assm+2][N_T+2];         // non-use
double duG_topdown[N_assm+2][N_T+2];    // non-use
double duG_G[N_assm+2][N_T+2];          // non-use

double drPY1[N_assm+2][N_T+2];          // for fraction of opened intrasyn. AMPA receptors 
double drSB1[N_assm+2][N_T+2];          // non-use
double dsF[N_assm+2][N_T+2];;           // non-use
double drLB1[N_assm+2][N_T+2];          
double AMPA_c1[N_assm+2][N_T+2];        
double GABASB_c1[N_assm+2][N_T+2][PERIOD+2];       // non-use
double GABALB_c1[N_assm+2][N_T+2][PERIOD+2];      

double w_rec_1[N_assm+2][N_T+2][N_assm+2][N_T+2];  // P-to-P synaptic connection weight in N_S
double w_fed_1[N_assm+2][N_T+2];                   // non-use
double w_lat_1[N_assm+2][N_T+2][N_T+2];            // B-t-P synaptic connection weight
double w_v1Ia_v2[N_assm+2][N_T+2][N_assm+2][N_T+2];// non-use
double w_v1Ib_v2[N_assm+2][N_T+2][N_assm+2][N_T+2];// non-use
double w_v1P_v2[N_assm+2][N_T+2][N_assm+2][N_T+2]; // topdown P(N_M)-to-P(N_S) synaptic connection weight
double w_v1G_v2[N_assm+2][N_T+2][N_assm+2][N_T+2]; // non-use
double wSB_PY1[N_assm+2][N_T+2];                   // non-use
double wLB_PY1[N_assm+2][N_assm+2][N_T+2];         // P-t-B synaptic connection weight
double wG_P[N_assm+2];                             // non-use
double wG_Ia[N_assm+2];                            // non-use
double difuPY1(int,int);                           // calculation in membrane pot.   
double difuSB1(int,int);                           // non-use  
double difuLB1(int,int);     
double difuG(int,int);                             // non-use 
double difrPY1(int,int);                           // calculation in membrane pot. in fraction of open channes
double difrSB1(int,int);                           // non-use
double difsF(int,int);                             // non-use
double difrLB1(int,int);   

// N_M
double uPY2[N_assm+2][N_T+2];   
double uSB2[N_assm+2][N_T+2];  
double uLB2[N_assm+2][N_T+2];  
double vPY2[N_assm+2][N_T+2][PERIOD+2];  
double vSB2[N_assm+2][N_T+2][PERIOD+2];  
double vLB2[N_assm+2][N_T+2][PERIOD+2];  
double rPY2[N_assm+2][N_T+2];  
double rPY2d[N_assm+2][N_T+2][PERIOD+2];  
double rSB2[N_assm+2][N_T+2];  
double rLB2[N_assm+2][N_T+2]; 
double GlutPY_c2[N_assm+2][N_T+2];    
double GlutSB_c2[N_assm+2][N_T+2];    
double GlutLB_c2[N_assm+2][N_T+2];    
double duPY_leak2[N_assm+2][N_T+2];  
double duPY_rec_2[N_assm+2][N_T+2];  
double duPY_fed_2[N_assm+2][N_T+2];  
double duPY_lat_2[N_assm+2][N_T+2];  
double duPY_ext_2[N_assm+2][N_T+2];  
double I_leak2; 
double I_rec_2; 
double I_fed_2; 
double I_lat_2; 
double I_ext_2; 

double duSB_leak2[N_assm+2][N_T+2]; 
double duSB_2[N_assm+2][N_T+2];
double duSB_ext_2;                   
double I_leakSB2; 
double I_fed_SB2; 
double I_ext_SB2; 
double duLB_leak2[N_assm+2][N_T+2];      
double duLB2[N_assm+2][N_T+2];
double duLB2_bottom_up[N_assm+2][N_T+2]; 
double duLB_ext2;                        
double I_leakLB2; 
double I_fed_LB2; 
double I_ext_LB2; 
double drPY2[N_assm+2][N_T+2]; 
double drSB2[N_assm+2][N_T+2];  
double drLB2[N_assm+2][N_T+2]; 
double AMPA_c2[N_assm+2][N_T+2];  
double GABASB_c2[N_assm+2][N_T+2][PERIOD+2]; 
double GABALB_c2[N_assm+2][N_T+2][PERIOD+2];   
double w_rec_2[N_assm+2][N_T+2][N_assm+2][N_T+2];            
double w_fed_2[N_assm+2][N_T+2];                  
double w_lat_2[N_assm+2][N_T+2][N_T+2]; 
double w_v2_v1[N_assm+2][N_T+2][N_assm+2][N_T+2] ;
double wSB_PY2[N_assm+2][N_T+2];                   
double wLB_PY2[N_assm+2][N_assm+2][N_T+2];         
double difuPY2(int,int);   
double difuSB2(int,int);    
double difuLB2(int,int);   
double difrPY2(int,int);  
double difrSB2(int,int);   
double difrLB2(int,int);   
double sigmoidPY(double); 
double sigmoidPY2(double); 
double sigmoidSB(double); 
double sigmoidSB2(double); 
double sigmoidLB(double); 
double sigmoidLB2(double); 

double OFSWCA=0.01;                   // ofset for output file within cell assembly
double OFSBCA=0.12;                   // ofset for output file between cell assembly
int    xi[N_assm+2][N_assm+2][N_T+2];
double rand01(long int *);            // radom number [0, 1.0]
void display(void);                   // display
void init(void);                      // initialization  
long int SEEDMP=1000;  

FILE *gaba0_01,*gaba0_11,*gaba0_21,*gaba0_31,*gaba0_41,*gaba0_51,*gaba0_61,*gaba0_71,*gaba0_81,*gaba0_91;
FILE *gaba1_01,*gaba1_11,*gaba1_21,*gaba1_31,*gaba1_41,*gaba1_51,*gaba1_61,*gaba1_71,*gaba1_81,*gaba1_91;
FILE *gaba2_01,*gaba2_11,*gaba2_21,*gaba2_31,*gaba2_41,*gaba2_51,*gaba2_61,*gaba2_71,*gaba2_81,*gaba2_91;
FILE *gaba3_01,*gaba3_11,*gaba3_21,*gaba3_31,*gaba3_41,*gaba3_51,*gaba3_61,*gaba3_71,*gaba3_81,*gaba3_91;
FILE *gaba4_01,*gaba4_11,*gaba4_21,*gaba4_31,*gaba4_41,*gaba4_51,*gaba4_61,*gaba4_71,*gaba4_81,*gaba4_91;
FILE *gaba5_01,*gaba5_11,*gaba5_21,*gaba5_31,*gaba5_41,*gaba5_51,*gaba5_61,*gaba5_71,*gaba5_81,*gaba5_91;
FILE *gaba6_01,*gaba6_11,*gaba6_21,*gaba6_31,*gaba6_41,*gaba6_51,*gaba6_61,*gaba6_71,*gaba6_81,*gaba6_91;
FILE *gaba7_01,*gaba7_11,*gaba7_21,*gaba7_31,*gaba7_41,*gaba7_51,*gaba7_61,*gaba7_71,*gaba7_81,*gaba7_91;

// N_S
FILE *uPY0_01,*uPY0_11,*uPY0_21,*uPY0_31,*uPY0_41,*uPY0_51,*uPY0_61,*uPY0_71,*uPY0_81,*uPY0_91;
FILE *uPY1_01,*uPY1_11,*uPY1_21,*uPY1_31,*uPY1_41,*uPY1_51,*uPY1_61,*uPY1_71,*uPY1_81,*uPY1_91;
FILE *uPY2_01,*uPY2_11,*uPY2_21,*uPY2_31,*uPY2_41,*uPY2_51,*uPY2_61,*uPY2_71,*uPY2_81,*uPY2_91;
FILE *uPY3_01,*uPY3_11,*uPY3_21,*uPY3_31,*uPY3_41,*uPY3_51,*uPY3_61,*uPY3_71,*uPY3_81,*uPY3_91;
FILE *uPY4_01,*uPY4_11,*uPY4_21,*uPY4_31,*uPY4_41,*uPY4_51,*uPY4_61,*uPY4_71,*uPY4_81,*uPY4_91;
FILE *uPY5_01,*uPY5_11,*uPY5_21,*uPY5_31,*uPY5_41,*uPY5_51,*uPY5_61,*uPY5_71,*uPY5_81,*uPY5_91;
FILE *uPY6_01,*uPY6_11,*uPY6_21,*uPY6_31,*uPY6_41,*uPY6_51,*uPY6_61,*uPY6_71,*uPY6_81,*uPY6_91;
FILE *uPY7_01,*uPY7_11,*uPY7_21,*uPY7_31,*uPY7_41,*uPY7_51,*uPY7_61,*uPY7_71,*uPY7_81,*uPY7_91;

FILE *uLB0_01,*uLB1_01,*uLB2_01,*uLB3_01,*uLB4_01,*uLB5_01,*uLB6_01,*uLB7_01;

FILE *vLB0_01,*vLB0_11,*vLB0_21,*vLB0_31,*vLB0_41,*vLB0_51,*vLB0_61,*vLB0_71,*vLB0_81,*vLB0_91;
FILE *vLB1_01,*vLB1_11,*vLB1_21,*vLB1_31,*vLB1_41,*vLB1_51,*vLB1_61,*vLB1_71,*vLB1_81,*vLB1_91;
FILE *vLB2_01,*vLB2_11,*vLB2_21,*vLB2_31,*vLB2_41,*vLB2_51,*vLB2_61,*vLB2_71,*vLB2_81,*vLB2_91;
FILE *vLB3_01,*vLB3_11,*vLB3_21,*vLB3_31,*vLB3_41,*vLB3_51,*vLB3_61,*vLB3_71,*vLB3_81,*vLB3_91;
FILE *vLB4_01,*vLB4_11,*vLB4_21,*vLB4_31,*vLB4_41,*vLB4_51,*vLB4_61,*vLB4_71,*vLB4_81,*vLB4_91;
FILE *vLB5_01,*vLB5_11,*vLB5_21,*vLB5_31,*vLB5_41,*vLB5_51,*vLB5_61,*vLB5_71,*vLB5_81,*vLB5_91;
FILE *vLB6_01,*vLB6_11,*vLB6_21,*vLB6_31,*vLB6_41,*vLB6_51,*vLB6_61,*vLB6_71,*vLB6_81,*vLB6_91;
FILE *vLB7_01,*vLB7_11,*vLB7_21,*vLB7_31,*vLB7_41,*vLB7_51,*vLB7_61,*vLB7_71,*vLB7_81,*vLB7_91;

// N_M
FILE *uPY0_02,*uPY0_12,*uPY0_22,*uPY0_32,*uPY0_42,*uPY0_52,*uPY0_62,*uPY0_72,*uPY0_82,*uPY0_92;
FILE *uPY1_02,*uPY1_12,*uPY1_22,*uPY1_32,*uPY1_42,*uPY1_52,*uPY1_62,*uPY1_72,*uPY1_82,*uPY1_92;
FILE *uPY2_02,*uPY2_12,*uPY2_22,*uPY2_32,*uPY2_42,*uPY2_52,*uPY2_62,*uPY2_72,*uPY2_82,*uPY2_92;
FILE *uPY3_02,*uPY3_12,*uPY3_22,*uPY3_32,*uPY3_42,*uPY3_52,*uPY3_62,*uPY3_72,*uPY3_82,*uPY3_92;
FILE *uPY4_02,*uPY4_12,*uPY4_22,*uPY4_32,*uPY4_42,*uPY4_52,*uPY4_62,*uPY4_72,*uPY4_82,*uPY4_92;
FILE *uPY5_02,*uPY5_12,*uPY5_22,*uPY5_32,*uPY5_42,*uPY5_52,*uPY5_62,*uPY5_72,*uPY5_82,*uPY5_92;
FILE *uPY6_02,*uPY6_12,*uPY6_22,*uPY6_32,*uPY6_42,*uPY6_52,*uPY6_62,*uPY6_72,*uPY6_82,*uPY6_92;
FILE *uPY7_02,*uPY7_12,*uPY7_22,*uPY7_32,*uPY7_42,*uPY7_52,*uPY7_62,*uPY7_72,*uPY7_82,*uPY7_92;
FILE *leak2,*rec2,*fed2,*lat2,*ext2;  

FILE *uLB0_02,*uLB1_02,*uLB2_02,*uLB3_02,*uLB4_02,*uLB5_02,*uLB6_02,*uLB7_02;

FILE *vLB0_02,*vLB0_12,*vLB0_22,*vLB0_32,*vLB0_42,*vLB0_52,*vLB0_62,*vLB0_72,*vLB0_82,*vLB0_92;
FILE *vLB1_02,*vLB1_12,*vLB1_22,*vLB1_32,*vLB1_42,*vLB1_52,*vLB1_62,*vLB1_72,*vLB1_82,*vLB1_92;
FILE *vLB2_02,*vLB2_12,*vLB2_22,*vLB2_32,*vLB2_42,*vLB2_52,*vLB2_62,*vLB2_72,*vLB2_82,*vLB2_92;
FILE *vLB3_02,*vLB3_12,*vLB3_22,*vLB3_32,*vLB3_42,*vLB3_52,*vLB3_62,*vLB3_72,*vLB3_82,*vLB3_92;
FILE *vLB4_02,*vLB4_12,*vLB4_22,*vLB4_32,*vLB4_42,*vLB4_52,*vLB4_62,*vLB4_72,*vLB4_82,*vLB4_92;
FILE *vLB5_02,*vLB5_12,*vLB5_22,*vLB5_32,*vLB5_42,*vLB5_52,*vLB5_62,*vLB5_72,*vLB5_82,*vLB5_92;
FILE *vLB6_02,*vLB6_12,*vLB6_22,*vLB6_32,*vLB6_42,*vLB6_52,*vLB6_62,*vLB6_72,*vLB6_82,*vLB6_92;
FILE *vLB7_02,*vLB7_12,*vLB7_22,*vLB7_32,*vLB7_42,*vLB7_52,*vLB7_62,*vLB7_72,*vLB7_82,*vLB7_92;

void main(void){

	int theta,i,ii;
	double sigPYT,sigSBT,sigLBT;
	int ActPY1_ON[N_assm+2][N_T+2],ActSB1_ON[N_assm+2][N_T+2],ActLB1_ON[N_assm+2][N_T+2],duration;
	int ActPY2_ON[N_assm+2][N_T+2],ActSB2_ON[N_assm+2][N_T+2],ActLB2_ON[N_assm+2][N_T+2];
	int ActPY1_OF[N_assm+2][N_T+2],ActSB1_OF[N_assm+2][N_T+2],ActLB1_OF[N_assm+2][N_T+2];
	int ActPY2_OF[N_assm+2][N_T+2],ActSB2_OF[N_assm+2][N_T+2],ActLB2_OF[N_assm+2][N_T+2];

	init(); 

	gaba0_01=fopen("gaba0_01.dat","w");
	gaba0_11=fopen("gaba0_11.dat","w");
	gaba0_21=fopen("gaba0_21.dat","w");
	gaba0_31=fopen("gaba0_31.dat","w");
	gaba0_41=fopen("gaba0_41.dat","w");
	gaba0_51=fopen("gaba0_51.dat","w");
	gaba0_61=fopen("gaba0_61.dat","w");
	gaba0_71=fopen("gaba0_71.dat","w");
	gaba0_81=fopen("gaba0_81.dat","w");
	gaba0_91=fopen("gaba0_91.dat","w");
	gaba1_01=fopen("gaba1_01.dat","w");
	gaba1_11=fopen("gaba1_11.dat","w");
	gaba1_21=fopen("gaba1_21.dat","w");
	gaba1_31=fopen("gaba1_31.dat","w");
	gaba1_41=fopen("gaba1_41.dat","w");
	gaba1_51=fopen("gaba1_51.dat","w");
	gaba1_61=fopen("gaba1_61.dat","w");
	gaba1_71=fopen("gaba1_71.dat","w");
	gaba1_81=fopen("gaba1_81.dat","w");
	gaba1_91=fopen("gaba1_91.dat","w");
	gaba2_01=fopen("gaba2_01.dat","w");
	gaba2_11=fopen("gaba2_11.dat","w");
	gaba2_21=fopen("gaba2_21.dat","w");
	gaba2_31=fopen("gaba2_31.dat","w");
	gaba2_41=fopen("gaba2_41.dat","w");
	gaba2_51=fopen("gaba2_51.dat","w");
	gaba2_61=fopen("gaba2_61.dat","w");
	gaba2_71=fopen("gaba2_71.dat","w");
	gaba2_81=fopen("gaba2_81.dat","w");
	gaba2_91=fopen("gaba2_91.dat","w");
	gaba3_01=fopen("gaba3_01.dat","w");
	gaba3_11=fopen("gaba3_11.dat","w");
	gaba3_21=fopen("gaba3_21.dat","w");
	gaba3_31=fopen("gaba3_31.dat","w");
	gaba3_41=fopen("gaba3_41.dat","w");
	gaba3_51=fopen("gaba3_51.dat","w");
	gaba3_61=fopen("gaba3_61.dat","w");
	gaba3_71=fopen("gaba3_71.dat","w");
	gaba3_81=fopen("gaba3_81.dat","w");
	gaba3_91=fopen("gaba3_91.dat","w");
	gaba4_01=fopen("gaba4_01.dat","w");
	gaba4_11=fopen("gaba4_11.dat","w");
	gaba4_21=fopen("gaba4_21.dat","w");
	gaba4_31=fopen("gaba4_31.dat","w");
	gaba4_41=fopen("gaba4_41.dat","w");
	gaba4_51=fopen("gaba4_51.dat","w");
	gaba4_61=fopen("gaba4_61.dat","w");
	gaba4_71=fopen("gaba4_71.dat","w");
	gaba4_81=fopen("gaba4_81.dat","w");
	gaba4_91=fopen("gaba4_91.dat","w");
	gaba5_01=fopen("gaba5_01.dat","w");
	gaba5_11=fopen("gaba5_11.dat","w");
	gaba5_21=fopen("gaba5_21.dat","w");
	gaba5_31=fopen("gaba5_31.dat","w");
	gaba5_41=fopen("gaba5_41.dat","w");
	gaba5_51=fopen("gaba5_51.dat","w");
	gaba5_61=fopen("gaba5_61.dat","w");
	gaba5_71=fopen("gaba5_71.dat","w");
	gaba5_81=fopen("gaba5_81.dat","w");
	gaba5_91=fopen("gaba5_91.dat","w");
	gaba6_01=fopen("gaba6_01.dat","w");
	gaba6_11=fopen("gaba6_11.dat","w");
	gaba6_21=fopen("gaba6_21.dat","w");
	gaba6_31=fopen("gaba6_31.dat","w");
	gaba6_41=fopen("gaba6_41.dat","w");
	gaba6_51=fopen("gaba6_51.dat","w");
	gaba6_61=fopen("gaba6_61.dat","w");
	gaba6_71=fopen("gaba6_71.dat","w");
	gaba6_81=fopen("gaba6_81.dat","w");
	gaba6_91=fopen("gaba6_91.dat","w");
	gaba7_01=fopen("gaba7_01.dat","w");
	gaba7_11=fopen("gaba7_11.dat","w");
	gaba7_21=fopen("gaba7_21.dat","w");
	gaba7_31=fopen("gaba7_31.dat","w");
	gaba7_41=fopen("gaba7_41.dat","w");
	gaba7_51=fopen("gaba7_51.dat","w");
	gaba7_61=fopen("gaba7_61.dat","w");
	gaba7_71=fopen("gaba7_71.dat","w");
	gaba7_81=fopen("gaba7_81.dat","w");
	gaba7_91=fopen("gaba7_91.dat","w");
    
	uPY0_01=fopen("uPY0_01.dat","w");
	uPY0_11=fopen("uPY0_11.dat","w");
	uPY0_21=fopen("uPY0_21.dat","w");
	uPY0_31=fopen("uPY0_31.dat","w");
	uPY0_41=fopen("uPY0_41.dat","w");
	uPY0_51=fopen("uPY0_51.dat","w");
	uPY0_61=fopen("uPY0_61.dat","w");
	uPY0_71=fopen("uPY0_71.dat","w");
	uPY0_81=fopen("uPY0_81.dat","w");
	uPY0_91=fopen("uPY0_91.dat","w");
	uPY1_01=fopen("uPY1_01.dat","w");
	uPY1_11=fopen("uPY1_11.dat","w");
	uPY1_21=fopen("uPY1_21.dat","w");
	uPY1_31=fopen("uPY1_31.dat","w");
	uPY1_41=fopen("uPY1_41.dat","w");
	uPY1_51=fopen("uPY1_51.dat","w");
	uPY1_61=fopen("uPY1_61.dat","w");
	uPY1_71=fopen("uPY1_71.dat","w");
	uPY1_81=fopen("uPY1_81.dat","w");
	uPY1_91=fopen("uPY1_91.dat","w");
	uPY2_01=fopen("uPY2_01.dat","w");
	uPY2_11=fopen("uPY2_11.dat","w");
	uPY2_21=fopen("uPY2_21.dat","w");
	uPY2_31=fopen("uPY2_31.dat","w");
	uPY2_41=fopen("uPY2_41.dat","w");
	uPY2_51=fopen("uPY2_51.dat","w");
	uPY2_61=fopen("uPY2_61.dat","w");
	uPY2_71=fopen("uPY2_71.dat","w");
	uPY2_81=fopen("uPY2_81.dat","w");
	uPY2_91=fopen("uPY2_91.dat","w");
	uPY3_01=fopen("uPY3_01.dat","w");
	uPY3_11=fopen("uPY3_11.dat","w");
	uPY3_21=fopen("uPY3_21.dat","w");
	uPY3_31=fopen("uPY3_31.dat","w");
	uPY3_41=fopen("uPY3_41.dat","w");
	uPY3_51=fopen("uPY3_51.dat","w");
	uPY3_61=fopen("uPY3_61.dat","w");
	uPY3_71=fopen("uPY3_71.dat","w");
	uPY3_81=fopen("uPY3_81.dat","w");
	uPY3_91=fopen("uPY3_91.dat","w");
	uPY4_01=fopen("uPY4_01.dat","w");
	uPY4_11=fopen("uPY4_11.dat","w");
	uPY4_21=fopen("uPY4_21.dat","w");
	uPY4_31=fopen("uPY4_31.dat","w");
	uPY4_41=fopen("uPY4_41.dat","w");
	uPY4_51=fopen("uPY4_51.dat","w");
	uPY4_61=fopen("uPY4_61.dat","w");
	uPY4_71=fopen("uPY4_71.dat","w");
	uPY4_81=fopen("uPY4_81.dat","w");
	uPY4_91=fopen("uPY4_91.dat","w");
	uPY5_01=fopen("uPY5_01.dat","w");
	uPY5_11=fopen("uPY5_11.dat","w");
	uPY5_21=fopen("uPY5_21.dat","w");
	uPY5_31=fopen("uPY5_31.dat","w");
	uPY5_41=fopen("uPY5_41.dat","w");
	uPY5_51=fopen("uPY5_51.dat","w");
	uPY5_61=fopen("uPY5_61.dat","w");
	uPY5_71=fopen("uPY5_71.dat","w");
	uPY5_81=fopen("uPY5_81.dat","w");
	uPY5_91=fopen("uPY5_91.dat","w");
	uPY6_01=fopen("uPY6_01.dat","w");
	uPY6_11=fopen("uPY6_11.dat","w");
	uPY6_21=fopen("uPY6_21.dat","w");
	uPY6_31=fopen("uPY6_31.dat","w");
	uPY6_41=fopen("uPY6_41.dat","w");
	uPY6_51=fopen("uPY6_51.dat","w");
	uPY6_61=fopen("uPY6_61.dat","w");
	uPY6_71=fopen("uPY6_71.dat","w");
	uPY6_81=fopen("uPY6_81.dat","w");
	uPY6_91=fopen("uPY6_91.dat","w");
	uPY7_01=fopen("uPY7_01.dat","w");
	uPY7_11=fopen("uPY7_11.dat","w");
	uPY7_21=fopen("uPY7_21.dat","w");
	uPY7_31=fopen("uPY7_31.dat","w");
	uPY7_41=fopen("uPY7_41.dat","w");
	uPY7_51=fopen("uPY7_51.dat","w");
	uPY7_61=fopen("uPY7_61.dat","w");
	uPY7_71=fopen("uPY7_71.dat","w");
	uPY7_81=fopen("uPY7_81.dat","w");
	uPY7_91=fopen("uPY7_91.dat","w");

	uLB0_01=fopen("uLB0_01.dat","w");
	uLB1_01=fopen("uLB1_01.dat","w");
	uLB2_01=fopen("uLB2_01.dat","w");
	uLB3_01=fopen("uLB3_01.dat","w");
	uLB4_01=fopen("uLB4_01.dat","w");
	uLB5_01=fopen("uLB5_01.dat","w");
	uLB6_01=fopen("uLB6_01.dat","w");
	uLB7_01=fopen("uLB7_01.dat","w");

	uPY0_02=fopen("uPY0_02.dat","w");
	uPY0_12=fopen("uPY0_12.dat","w");
	uPY0_22=fopen("uPY0_22.dat","w");
	uPY0_32=fopen("uPY0_32.dat","w");
	uPY0_42=fopen("uPY0_42.dat","w");
	uPY0_52=fopen("uPY0_52.dat","w");
	uPY0_62=fopen("uPY0_62.dat","w");
	uPY0_72=fopen("uPY0_72.dat","w");
	uPY0_82=fopen("uPY0_82.dat","w");
	uPY0_92=fopen("uPY0_92.dat","w");
	uPY1_02=fopen("uPY1_02.dat","w");
	uPY1_12=fopen("uPY1_12.dat","w");
	uPY1_22=fopen("uPY1_22.dat","w");
	uPY1_32=fopen("uPY1_32.dat","w");
	uPY1_42=fopen("uPY1_42.dat","w");
	uPY1_52=fopen("uPY1_52.dat","w");
	uPY1_62=fopen("uPY1_62.dat","w");
	uPY1_72=fopen("uPY1_72.dat","w");
	uPY1_82=fopen("uPY1_82.dat","w");
	uPY1_92=fopen("uPY1_92.dat","w");
	uPY2_02=fopen("uPY2_02.dat","w");
	uPY2_12=fopen("uPY2_12.dat","w");
	uPY2_22=fopen("uPY2_22.dat","w");
	uPY2_32=fopen("uPY2_32.dat","w");
	uPY2_42=fopen("uPY2_42.dat","w");
	uPY2_52=fopen("uPY2_52.dat","w");
	uPY2_62=fopen("uPY2_62.dat","w");
	uPY2_72=fopen("uPY2_72.dat","w");
	uPY2_82=fopen("uPY2_82.dat","w");
	uPY2_92=fopen("uPY2_92.dat","w");
	uPY3_02=fopen("uPY3_02.dat","w");
	uPY3_12=fopen("uPY3_12.dat","w");
	uPY3_22=fopen("uPY3_22.dat","w");
	uPY3_32=fopen("uPY3_32.dat","w");
	uPY3_42=fopen("uPY3_42.dat","w");
	uPY3_52=fopen("uPY3_52.dat","w");
	uPY3_62=fopen("uPY3_62.dat","w");
	uPY3_72=fopen("uPY3_72.dat","w");
	uPY3_82=fopen("uPY3_82.dat","w");
	uPY3_92=fopen("uPY3_92.dat","w");
	
	uPY4_02=fopen("uPY4_02.dat","w");
	uPY4_12=fopen("uPY4_12.dat","w");
	uPY4_22=fopen("uPY4_22.dat","w");
	uPY4_32=fopen("uPY4_32.dat","w");
	uPY4_42=fopen("uPY4_42.dat","w");
	uPY4_52=fopen("uPY4_52.dat","w");
	uPY4_62=fopen("uPY4_62.dat","w");
	uPY4_72=fopen("uPY4_72.dat","w");
	uPY4_82=fopen("uPY4_82.dat","w");
	uPY4_92=fopen("uPY4_92.dat","w");
	uPY5_02=fopen("uPY5_02.dat","w");
	uPY5_12=fopen("uPY5_12.dat","w");
	uPY5_22=fopen("uPY5_22.dat","w");
	uPY5_32=fopen("uPY5_32.dat","w");
	uPY5_42=fopen("uPY5_42.dat","w");
	uPY5_52=fopen("uPY5_52.dat","w");
	uPY5_62=fopen("uPY5_62.dat","w");
	uPY5_72=fopen("uPY5_72.dat","w");
	uPY5_82=fopen("uPY5_82.dat","w");
	uPY5_92=fopen("uPY5_92.dat","w");
	uPY6_02=fopen("uPY6_02.dat","w");
	uPY6_12=fopen("uPY6_12.dat","w");
	uPY6_22=fopen("uPY6_22.dat","w");
	uPY6_32=fopen("uPY6_32.dat","w");
	uPY6_42=fopen("uPY6_42.dat","w");
	uPY6_52=fopen("uPY6_52.dat","w");
	uPY6_62=fopen("uPY6_62.dat","w");
	uPY6_72=fopen("uPY6_72.dat","w");
	uPY6_82=fopen("uPY6_82.dat","w");
	uPY6_92=fopen("uPY6_92.dat","w");
	uPY7_02=fopen("uPY7_02.dat","w");
	uPY7_12=fopen("uPY7_12.dat","w");
	uPY7_22=fopen("uPY7_22.dat","w");
	uPY7_32=fopen("uPY7_32.dat","w");
	uPY7_42=fopen("uPY7_42.dat","w");
	uPY7_52=fopen("uPY7_52.dat","w");
	uPY7_62=fopen("uPY7_62.dat","w");
	uPY7_72=fopen("uPY7_72.dat","w");
	uPY7_82=fopen("uPY7_82.dat","w");
	uPY7_92=fopen("uPY7_92.dat","w");
	
	uLB0_02=fopen("uLB0_02.dat","w");
	uLB1_02=fopen("uLB1_02.dat","w");
	uLB2_02=fopen("uLB2_02.dat","w");
	uLB3_02=fopen("uLB3_02.dat","w");
	uLB4_02=fopen("uLB4_02.dat","w");
	uLB5_02=fopen("uLB5_02.dat","w");
	uLB6_02=fopen("uLB6_02.dat","w");
	uLB7_02=fopen("uLB7_02.dat","w");

	for (theta=0; theta<=N_assm; ++theta){ 
		for (ii=0; ii<=N_T; ++ii){
			ActPY1_ON[theta][ii] = 0; 
			ActSB1_ON[theta][ii] = 0; 
			ActLB1_ON[theta][ii] = 0; 
			ActPY1_OF[theta][ii] = 0; 
			ActSB1_OF[theta][ii] = 0; 
			ActLB1_OF[theta][ii] = 0; 
			uPY1[theta][ii]=UPYres;
			uSB1[theta][ii]=USBres;
			uLB1[theta][ii]=ULBres;
		}
	}
	for (theta=0; theta<=N_assm; ++theta){  
		for (ii=0; ii<=N_T; ++ii){
			ActPY2_ON[theta][ii] = 0; 
			ActSB2_ON[theta][ii] = 0; 
			ActLB2_ON[theta][ii] = 0; 
			ActPY2_OF[theta][ii] = 0; 
			ActSB2_OF[theta][ii] = 0; 
			ActLB2_OF[theta][ii] = 0; 
			uPY2[theta][ii]=UPYres;
			uSB2[theta][ii]=USBres;
			uLB2[theta][ii]=ULBres;
		}
	}
	duration  = (int)(0.001/DT);

	for (t=0; t<PERIOD; ++t){    
		display();
		for (theta=0; theta<=N_assm-1; ++theta){
			for (i=0; i<=N_T; ++i){	
				if (theta==COLUMN && i==NEURON){  
					I_leak1  = duPY_leak1[theta][i]*(cmPY/DT);
					I_rec_1  = duPY_rec_1[theta][i]*(cmPY/DT);
					I_fed_1  = duPY_fed_1[theta][i]*(cmPY/DT);
					I_lat_1  = duPY_lat_1[theta][i]*(cmPY/DT);
					I_topdown= duPY_topdown[theta][i]*(cmPY/DT);
					I_ext_1  = duPY_ext_1[theta][i]*(cmPY/DT);
					I_MG1    = duPY_MGB1[theta]*(cmPY/DT)*(int)(rand()/32768.0+inp_prob);
				}
				uPY1[theta][i] += difuPY1(theta,i); 
				sigPYT = sigmoidPY(uPY1[theta][i]);
				if (ActPY1_ON[theta][i]==0 && ActPY1_OF[theta][i]==0) vPY1[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigPYT);
				if (vPY1[theta][i][t]==1.0){
					uPY1[theta][i]=UPYact;
					GlutPY_c1[theta][i]=Glut_c;
					ActPY1_ON[theta][i]=1;
					ActPY1_OF[theta][i]=1;
				}else{
					GlutPY_c1[theta][i]=0.0;
				}
				if (ActPY1_ON[theta][i] != 0){
					++ActPY1_ON[theta][i];
					uPY1[theta][i]=UPYact;
					GlutPY_c1[theta][i]=Glut_c;
				}
				if (ActPY1_OF[theta][i] != 0){
					++ActPY1_OF[theta][i];
				}
				if (ActPY1_ON[theta][i] > duration){
					uPY1[theta][i]=UPYres;
					ActPY1_ON[theta][i] = 0;
				}
				if (ActPY1_OF[theta][i] > 10*duration){
					ActPY1_OF[theta][i] = 0;
				}
				rPY1[theta][i] += difrPY1(theta,i);   
				rPY1d[theta][i][t] = rPY1[theta][i]; 
				if (t%INTV == 0){
				}

				uSB1[theta][i] += difuSB1(theta,i); 
				sigSBT = sigmoidSB(uSB1[theta][i]);
				if (ActSB1_ON[theta][i]==0 && ActSB1_OF[theta][i]==0) vSB1[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigSBT);
				if (vSB1[theta][i][t]==1.0){
					uSB1[theta][i]=USBact;
					GABASB_c1[theta][i][t]=GABA_cS;
					ActSB1_ON[theta][i]=1;
					ActSB1_OF[theta][i]=1;
				}else{
					GABASB_c1[theta][i][t]=0.0;
				}
				if (ActSB1_ON[theta][i] != 0){
					++ActSB1_ON[theta][i];
					uSB1[theta][i]=USBact;
					GABASB_c1[theta][i][t]=GABA_cS;
				}
				if (ActSB1_OF[theta][i] != 0){
					++ActSB1_OF[theta][i];
				}
				if (ActSB1_ON[theta][i] > duration){
					uSB1[theta][i]=USBres;
					ActSB1_ON[theta][i] = 0;
				}
				if (ActSB1_OF[theta][i] > 10*duration){
					ActSB1_OF[theta][i] = 0;
				}	
				rSB1[theta][i] += difrSB1(theta,i); 
				sF[theta][i]   += difsF(theta,i);      
				if (t%INTV == 0){
				}

				uLB1[theta][i] += difuLB1(theta,i); 
				sigLBT = sigmoidLB(uLB1[theta][i]);
				if (ActLB1_ON[theta][i]==0 && ActLB1_OF[theta][i]==0) vLB1[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigLBT);
				if (vLB1[theta][i][t]==1.0){
					uLB1[theta][i]=ULBact;
					GABALB_c1[theta][i][t]=GABA_cL1;
					ActLB1_ON[theta][i]=1;
					ActLB1_OF[theta][i]=1;
				}else{
					GABALB_c1[theta][i][t]=0.0;
				}
				if (ActLB1_ON[theta][i] != 0){
					++ActLB1_ON[theta][i];
					uLB1[theta][i]=ULBact;
					GABALB_c1[theta][i][t]=GABA_cL1;
				}
				if (ActLB1_OF[theta][i] != 0){
					++ActLB1_OF[theta][i];
				}				
				if (ActLB1_ON[theta][i] > duration){
					uLB1[theta][i]=ULBres;
					ActLB1_ON[theta][i] = 0;
				}
				if (ActLB1_OF[theta][i] > 10*duration){
					ActLB1_OF[theta][i] = 0;
				}
				rLB1[theta][i] += difrLB1(theta,i); 
				if (t%INTV == 0){
				}

				uG[theta][i] += difuG(theta,i); 

				if (theta==COLUMN && i==NEURON){  
					I_leak2 = duPY_leak2[theta][i]*(cmPY/DT);
					I_rec_2 = duPY_rec_2[theta][i]*(cmPY/DT);
					I_fed_2 = duPY_fed_2[theta][i]*(cmPY/DT);
					I_lat_2 = duPY_lat_2[theta][i]*(cmPY/DT);
					I_ext_2 = duPY_ext_2[theta][i]*(cmPY/DT);
				}
				uPY2[theta][i] += difuPY2(theta,i); 
				sigPYT = sigmoidPY2(uPY2[theta][i]);
				if (ActPY2_ON[theta][i]==0 && ActPY2_OF[theta][i]==0) vPY2[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigPYT);
				if (vPY2[theta][i][t]==1.0){
					uPY2[theta][i]=UPYact;
					GlutPY_c2[theta][i]=Glut_c;
					ActPY2_ON[theta][i]=1;
					ActPY2_OF[theta][i]=1;
				}else{
					GlutPY_c2[theta][i]=0.0;
				}
				if (ActPY2_ON[theta][i] != 0){
					++ActPY2_ON[theta][i];
					uPY2[theta][i]=UPYact;
					GlutPY_c2[theta][i]=Glut_c;
				}
				if (ActPY2_OF[theta][i] != 0){
					++ActPY2_OF[theta][i];
				}
				if (ActPY2_ON[theta][i] > duration){
					uPY2[theta][i]=UPYres;
					ActPY2_ON[theta][i] = 0;
				}
				if (ActPY2_OF[theta][i] > 10*duration){
					ActPY2_OF[theta][i] = 0;
				}
				rPY2[theta][i] += difrPY2(theta,i); 
				rPY2d[theta][i][t] = rPY2[theta][i];

				uSB2[theta][i] += difuSB2(theta,i); 
				sigSBT = sigmoidSB2(uSB2[theta][i]);
				if (ActSB2_ON[theta][i]==0 && ActSB2_OF[theta][i]==0) vSB2[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigSBT);
				if (vSB2[theta][i][t]==1.0){
					uSB2[theta][i]=USBact;
					GABASB_c2[theta][i][t]=GABA_cS;
					ActSB2_ON[theta][i]=1;
					ActSB2_OF[theta][i]=1;
				}else{
					GABASB_c2[theta][i][t]=0.0;
				}
				if (ActSB2_ON[theta][i] != 0){
					++ActSB2_ON[theta][i];
					uSB2[theta][i]=USBact;
					GABASB_c2[theta][i][t]=GABA_cS;
				}
				if (ActSB2_OF[theta][i] != 0){
					++ActSB2_OF[theta][i];
				}
				if (ActSB2_ON[theta][i] > duration){
					uSB2[theta][i]=USBres;
					ActSB2_ON[theta][i] = 0;
				}
				if (ActSB2_OF[theta][i] > 10*duration){
					ActSB2_OF[theta][i] = 0;
				}	
				rSB2[theta][i] += difrSB2(theta,i); 
				if (t%INTV == 0){
				}

				uLB2[theta][i] += difuLB2(theta,i); 
				sigLBT = sigmoidLB2(uLB2[theta][i]);
				if (ActLB2_ON[theta][i]==0 && ActLB2_OF[theta][i]==0) vLB2[theta][i][t]=(double)(int)(rand01(&SEEDMP)+sigLBT);
				if (vLB2[theta][i][t]==1.0){
					uLB2[theta][i]=ULBact;
					GABALB_c2[theta][i][t]=GABA_cL2;
					ActLB2_ON[theta][i]=1;
					ActLB2_OF[theta][i]=1;
				}else{
					GABALB_c2[theta][i][t]=0.0;
				}
				if (ActLB2_ON[theta][i] != 0){
					++ActLB2_ON[theta][i];
					uLB2[theta][i]=ULBact;
					GABALB_c2[theta][i][t]=GABA_cL2;
				}
				if (ActLB2_OF[theta][i] != 0){
					++ActLB2_OF[theta][i];
				}				
				if (ActLB2_ON[theta][i] > duration){
					uLB2[theta][i]=ULBres;
					ActLB2_ON[theta][i] = 0;
				}
				if (ActLB2_OF[theta][i] > 10*duration){
					ActLB2_OF[theta][i] = 0;
				}
				rLB2[theta][i] += difrLB2(theta,i); 
				if (t%INTV == 0){
				}

				if (t>=OUT){  
					if (theta==0 && i==0){									
						fprintf(uPY0_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY0_02,"%f\n",uPY2[theta][i]); 
						fprintf(uLB0_01,"%f\n",uLB1[theta][i]); 
						fprintf(gaba0_01,"%e\n",GABA_ext[theta][i]);
					} 
					if (theta==0 && i==1){									
						fprintf(uPY0_11,"%f\n",uPY1[theta][i]); 
						fprintf(uPY0_12,"%f\n",uPY2[theta][i]); 
						fprintf(gaba0_11,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==0 && i==2){									
						fprintf(uPY0_21,"%f\n",uPY1[theta][i]); 
						fprintf(uPY0_22,"%f\n",uPY2[theta][i]); 
						fprintf(gaba0_21,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==0 && i==3){									
						fprintf(uPY0_31,"%f\n",uPY1[theta][i]); 
						fprintf(uPY0_32,"%f\n",uPY2[theta][i]); 
						fprintf(gaba0_31,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==0 && i==4){									
						fprintf(uPY0_41,"%f\n",uPY1[theta][i]); 
						fprintf(uPY0_42,"%f\n",uPY2[theta][i]); 
						fprintf(gaba0_41,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==0 && i==5){									
						fprintf(uPY0_51,"%f\n",uPY1[theta][i]); 
						fprintf(uPY0_52,"%f\n",uPY2[theta][i]); 
						fprintf(gaba0_51,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==0 && i==6){									
						fprintf(uPY0_61,"%f\n",uPY1[theta][i]); 
						fprintf(uPY0_62,"%f\n",uPY2[theta][i]); 
						fprintf(gaba0_61,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==0 && i==7){									
						fprintf(uPY0_71,"%f\n",uPY1[theta][i]); 
						fprintf(uPY0_72,"%f\n",uPY2[theta][i]); 
						fprintf(gaba0_71,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==0 && i==8){									
						fprintf(uPY0_81,"%f\n",uPY1[theta][i]); 
						fprintf(uPY0_82,"%f\n",uPY2[theta][i]); 
						fprintf(gaba0_81,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==0 && i==9){									
						fprintf(uPY0_91,"%f\n",uPY1[theta][i]); 
						fprintf(uPY0_92,"%f\n",uPY2[theta][i]); 
						fprintf(gaba0_91,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==1 && i==0){									
						fprintf(uPY1_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY1_02,"%f\n",uPY2[theta][i]); 
						fprintf(uLB1_01,"%f\n",uLB1[theta][i]); 
						fprintf(uLB1_02,"%f\n",uLB2[theta][i]); 
						fprintf(gaba1_01,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==1 && i==1){									
						fprintf(uPY1_11,"%f\n",uPY1[theta][i]); 
						fprintf(uPY1_12,"%f\n",uPY2[theta][i]); 
						fprintf(gaba1_11,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==1 && i==2){									
						fprintf(uPY1_21,"%f\n",uPY1[theta][i]); 
						fprintf(uPY1_22,"%f\n",uPY2[theta][i]); 
						fprintf(gaba1_21,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==1 && i==3){									
						fprintf(uPY1_31,"%f\n",uPY1[theta][i]); 
						fprintf(uPY1_32,"%f\n",uPY2[theta][i]); 
						fprintf(gaba1_31,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==1 && i==4){									
						fprintf(uPY1_41,"%f\n",uPY1[theta][i]); 
						fprintf(uPY1_42,"%f\n",uPY2[theta][i]); 
						fprintf(gaba1_41,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==1 && i==5){									
						fprintf(uPY1_51,"%f\n",uPY1[theta][i]); 
						fprintf(uPY1_52,"%f\n",uPY2[theta][i]); 
						fprintf(gaba1_51,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==1 && i==6){									
						fprintf(uPY1_61,"%f\n",uPY1[theta][i]); 
						fprintf(uPY1_62,"%f\n",uPY2[theta][i]); 
						fprintf(gaba1_61,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==1 && i==7){									
						fprintf(uPY1_71,"%f\n",uPY1[theta][i]); 
						fprintf(uPY1_72,"%f\n",uPY2[theta][i]); 
						fprintf(gaba1_71,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==1 && i==8){									
						fprintf(uPY1_81,"%f\n",uPY1[theta][i]); 
						fprintf(uPY1_82,"%f\n",uPY2[theta][i]); 
						fprintf(gaba1_81,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==1 && i==9){									
						fprintf(uPY1_91,"%f\n",uPY1[theta][i]); 
						fprintf(uPY1_92,"%f\n",uPY2[theta][i]); 
						fprintf(gaba1_91,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==2 && i==0){									
						fprintf(uPY2_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY2_02,"%f\n",uPY2[theta][i]); 
						fprintf(uLB2_01,"%f\n",uLB1[theta][i]); 
						fprintf(uLB2_02,"%f\n",uLB2[theta][i]); 
						fprintf(gaba2_01,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==2 && i==1){									
						fprintf(uPY2_11,"%f\n",uPY1[theta][i]); 
						fprintf(uPY2_12,"%f\n",uPY2[theta][i]); 
						fprintf(gaba2_11,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==2 && i==2){									
						fprintf(uPY2_21,"%f\n",uPY1[theta][i]); 
						fprintf(uPY2_22,"%f\n",uPY2[theta][i]); 
						fprintf(gaba2_21,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==2 && i==3){									
						fprintf(uPY2_31,"%f\n",uPY1[theta][i]); 
						fprintf(uPY2_32,"%f\n",uPY2[theta][i]); 
						fprintf(gaba2_31,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==2 && i==4){									
						fprintf(uPY2_41,"%f\n",uPY1[theta][i]); 
						fprintf(uPY2_42,"%f\n",uPY2[theta][i]); 
						fprintf(gaba2_41,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==2 && i==5){									
						fprintf(uPY2_51,"%f\n",uPY1[theta][i]); 
						fprintf(uPY2_52,"%f\n",uPY2[theta][i]); 
						fprintf(gaba2_51,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==2 && i==6){									
						fprintf(uPY2_61,"%f\n",uPY1[theta][i]); 
						fprintf(uPY2_62,"%f\n",uPY2[theta][i]); 
						fprintf(gaba2_61,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==2 && i==7){									
						fprintf(uPY2_71,"%f\n",uPY1[theta][i]); 
						fprintf(uPY2_72,"%f\n",uPY2[theta][i]); 
						fprintf(gaba2_71,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==2 && i==8){									
						fprintf(uPY2_81,"%f\n",uPY1[theta][i]); 
						fprintf(uPY2_82,"%f\n",uPY2[theta][i]); 
						fprintf(gaba2_81,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==2 && i==9){									
						fprintf(uPY2_91,"%f\n",uPY1[theta][i]); 
						fprintf(uPY2_92,"%f\n",uPY2[theta][i]); 
						fprintf(gaba2_91,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==3 && i==0){									
						fprintf(uPY3_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY3_02,"%f\n",uPY2[theta][i]); 
						fprintf(uLB3_01,"%f\n",uLB1[theta][i]); 
						fprintf(uLB3_02,"%f\n",uLB2[theta][i]); 
						fprintf(gaba3_01,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==3 && i==1){									
						fprintf(uPY3_11,"%f\n",uPY1[theta][i]); 
						fprintf(uPY3_12,"%f\n",uPY2[theta][i]); 
						fprintf(gaba3_11,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==3 && i==2){									
						fprintf(uPY3_21,"%f\n",uPY1[theta][i]); 
						fprintf(uPY3_22,"%f\n",uPY2[theta][i]); 
						fprintf(gaba3_21,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==3 && i==3){									
						fprintf(uPY3_31,"%f\n",uPY1[theta][i]); 
						fprintf(uPY3_32,"%f\n",uPY2[theta][i]); 
						fprintf(gaba3_31,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==3 && i==4){									
						fprintf(uPY3_41,"%f\n",uPY1[theta][i]); 
						fprintf(uPY3_42,"%f\n",uPY2[theta][i]); 
						fprintf(gaba3_41,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==3 && i==5){									
						fprintf(uPY3_51,"%f\n",uPY1[theta][i]); 
						fprintf(uPY3_52,"%f\n",uPY2[theta][i]); 
						fprintf(gaba3_51,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==3 && i==6){									
						fprintf(uPY3_61,"%f\n",uPY1[theta][i]); 
						fprintf(uPY3_62,"%f\n",uPY2[theta][i]); 
						fprintf(gaba3_61,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==3 && i==7){									
						fprintf(uPY3_71,"%f\n",uPY1[theta][i]); 
						fprintf(uPY3_72,"%f\n",uPY2[theta][i]); 
						fprintf(gaba3_71,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==3 && i==8){									
						fprintf(uPY3_81,"%f\n",uPY1[theta][i]); 
						fprintf(uPY3_82,"%f\n",uPY2[theta][i]); 
						fprintf(gaba3_81,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==3 && i==9){									
						fprintf(uPY3_91,"%f\n",uPY1[theta][i]); 
						fprintf(uPY3_92,"%f\n",uPY2[theta][i]); 
						fprintf(gaba3_91,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==4 && i==0){									
						fprintf(uPY4_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY4_02,"%f\n",uPY2[theta][i]); 
						fprintf(uLB4_01,"%f\n",uLB1[theta][i]); 
						fprintf(uLB4_02,"%f\n",uLB2[theta][i]); 
						fprintf(gaba4_01,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==4 && i==1){									
						fprintf(uPY4_11,"%f\n",uPY1[theta][i]); 
						fprintf(uPY4_12,"%f\n",uPY2[theta][i]); 
						fprintf(gaba4_11,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==4 && i==2){									
						fprintf(uPY4_21,"%f\n",uPY1[theta][i]); 
						fprintf(uPY4_22,"%f\n",uPY2[theta][i]); 
						fprintf(gaba4_21,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==4 && i==3){									
						fprintf(uPY4_31,"%f\n",uPY1[theta][i]); 
						fprintf(uPY4_32,"%f\n",uPY2[theta][i]); 
						fprintf(gaba4_31,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==4 && i==4){									
						fprintf(uPY4_41,"%f\n",uPY1[theta][i]); 
						fprintf(uPY4_42,"%f\n",uPY2[theta][i]); 
						fprintf(gaba4_41,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==4 && i==5){									
						fprintf(uPY4_51,"%f\n",uPY1[theta][i]); 
						fprintf(uPY4_52,"%f\n",uPY2[theta][i]); 
						fprintf(gaba4_51,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==4 && i==6){									
						fprintf(uPY4_61,"%f\n",uPY1[theta][i]); 
						fprintf(uPY4_62,"%f\n",uPY2[theta][i]); 
						fprintf(gaba4_61,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==4 && i==7){									
						fprintf(uPY4_71,"%f\n",uPY1[theta][i]); 
						fprintf(uPY4_72,"%f\n",uPY2[theta][i]); 
						fprintf(gaba4_71,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==4 && i==8){									
						fprintf(uPY4_81,"%f\n",uPY1[theta][i]); 
						fprintf(uPY4_82,"%f\n",uPY2[theta][i]); 
						fprintf(gaba4_81,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==4 && i==9){									
						fprintf(uPY4_91,"%f\n",uPY1[theta][i]); 
						fprintf(uPY4_92,"%f\n",uPY2[theta][i]); 
						fprintf(gaba4_91,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==5 && i==0){									
						fprintf(uPY5_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY5_02,"%f\n",uPY2[theta][i]); 
						fprintf(gaba5_01,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==5 && i==1){									
						fprintf(uPY5_11,"%f\n",uPY1[theta][i]); 
						fprintf(uPY5_12,"%f\n",uPY2[theta][i]); 
						fprintf(gaba5_11,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==5 && i==2){									
						fprintf(uPY5_21,"%f\n",uPY1[theta][i]); 
						fprintf(uPY5_22,"%f\n",uPY2[theta][i]); 
						fprintf(gaba5_21,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==5 && i==3){									
						fprintf(uPY5_31,"%f\n",uPY1[theta][i]); 
						fprintf(uPY5_32,"%f\n",uPY2[theta][i]); 
						fprintf(gaba5_31,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==5 && i==4){									
						fprintf(uPY5_41,"%f\n",uPY1[theta][i]); 
						fprintf(uPY5_42,"%f\n",uPY2[theta][i]); 
						fprintf(gaba5_41,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==5 && i==5){									
						fprintf(uPY5_51,"%f\n",uPY1[theta][i]); 
						fprintf(uPY5_52,"%f\n",uPY2[theta][i]); 
						fprintf(gaba5_51,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==5 && i==6){									
						fprintf(uPY5_61,"%f\n",uPY1[theta][i]); 
						fprintf(uPY5_62,"%f\n",uPY2[theta][i]); 
						fprintf(gaba5_61,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==5 && i==7){									
						fprintf(uPY5_71,"%f\n",uPY1[theta][i]); 
						fprintf(uPY5_72,"%f\n",uPY2[theta][i]); 
						fprintf(gaba5_71,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==5 && i==8){									
						fprintf(uPY5_81,"%f\n",uPY1[theta][i]); 
						fprintf(uPY5_82,"%f\n",uPY2[theta][i]); 
						fprintf(gaba5_81,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==5 && i==9){									
						fprintf(uPY5_91,"%f\n",uPY1[theta][i]); 
						fprintf(uPY5_92,"%f\n",uPY2[theta][i]); 
						fprintf(gaba5_91,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==6 && i==0){									
						fprintf(uPY6_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY6_02,"%f\n",uPY2[theta][i]); 
						fprintf(uLB6_01,"%f\n",uLB1[theta][i]); 
						fprintf(gaba6_01,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==6 && i==1){									
						fprintf(uPY6_11,"%f\n",uPY1[theta][i]); 
						fprintf(uPY6_12,"%f\n",uPY2[theta][i]); 
						fprintf(gaba6_11,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==6 && i==2){									
						fprintf(uPY6_21,"%f\n",uPY1[theta][i]); 
						fprintf(uPY6_22,"%f\n",uPY2[theta][i]); 
						fprintf(gaba6_21,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==6 && i==3){									
						fprintf(uPY6_31,"%f\n",uPY1[theta][i]); 
						fprintf(uPY6_32,"%f\n",uPY2[theta][i]); 
						fprintf(gaba6_31,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==6 && i==4){									
						fprintf(uPY6_41,"%f\n",uPY1[theta][i]); 
						fprintf(uPY6_42,"%f\n",uPY2[theta][i]); 
						fprintf(gaba6_41,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==6 && i==5){									
						fprintf(uPY6_51,"%f\n",uPY1[theta][i]); 
						fprintf(uPY6_52,"%f\n",uPY2[theta][i]); 
						fprintf(gaba6_51,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==6 && i==6){									
						fprintf(uPY6_61,"%f\n",uPY1[theta][i]); 
						fprintf(uPY6_62,"%f\n",uPY2[theta][i]); 
						fprintf(gaba6_61,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==6 && i==7){									
						fprintf(uPY6_71,"%f\n",uPY1[theta][i]); 
						fprintf(uPY6_72,"%f\n",uPY2[theta][i]); 
						fprintf(gaba6_71,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==6 && i==8){									
						fprintf(uPY6_81,"%f\n",uPY1[theta][i]); 
						fprintf(uPY6_82,"%f\n",uPY2[theta][i]); 
						fprintf(gaba6_81,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==6 && i==9){									
						fprintf(uPY6_91,"%f\n",uPY1[theta][i]); 
						fprintf(uPY6_92,"%f\n",uPY2[theta][i]); 
						fprintf(gaba6_91,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==7 && i==0){									
						fprintf(uPY7_01,"%f\n",uPY1[theta][i]); 
						fprintf(uPY7_02,"%f\n",uPY2[theta][i]); 
						fprintf(uLB7_01,"%f\n",uLB1[theta][i]); 
						fprintf(gaba7_01,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==7 && i==1){									
						fprintf(uPY7_11,"%f\n",uPY1[theta][i]); 
						fprintf(uPY7_12,"%f\n",uPY2[theta][i]); 
						fprintf(gaba7_11,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==7 && i==2){									
						fprintf(uPY7_21,"%f\n",uPY1[theta][i]); 
						fprintf(uPY7_22,"%f\n",uPY2[theta][i]); 
						fprintf(gaba7_21,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==7 && i==3){									
						fprintf(uPY7_31,"%f\n",uPY1[theta][i]); 
						fprintf(uPY7_32,"%f\n",uPY2[theta][i]); 
						fprintf(gaba7_31,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==7 && i==4){									
						fprintf(uPY7_41,"%f\n",uPY1[theta][i]); 
						fprintf(uPY7_42,"%f\n",uPY2[theta][i]); 
						fprintf(gaba7_41,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==7 && i==5){									
						fprintf(uPY7_51,"%f\n",uPY1[theta][i]); 
						fprintf(uPY7_52,"%f\n",uPY2[theta][i]); 
						fprintf(gaba7_51,"%e\n",GABA_ext[theta][i]);
				}
					if (theta==7 && i==6){									
						fprintf(uPY7_61,"%f\n",uPY1[theta][i]); 
						fprintf(uPY7_62,"%f\n",uPY2[theta][i]); 
						fprintf(gaba7_61,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==7 && i==7){									
						fprintf(uPY7_71,"%f\n",uPY1[theta][i]); 
						fprintf(uPY7_72,"%f\n",uPY2[theta][i]); 
						fprintf(gaba7_71,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==7 && i==8){									
						fprintf(uPY7_81,"%f\n",uPY1[theta][i]); 
						fprintf(uPY7_82,"%f\n",uPY2[theta][i]); 
						fprintf(gaba7_81,"%e\n",GABA_ext[theta][i]);
					}
					if (theta==7 && i==9){									
						fprintf(uPY7_91,"%f\n",uPY1[theta][i]); 
						fprintf(uPY7_92,"%f\n",uPY2[theta][i]); 
						fprintf(gaba7_91,"%e\n",GABA_ext[theta][i]);
					}
				} 
			} 
		} 

		dfsGABA_ext(); 
		for (theta=0; theta<=N_assm; ++theta){
			for (i=0; i<=N_T; ++i){	
			rEXT[theta][i] += difrEXT(theta,i); 
			}
		}
		for (theta=0; theta<=N_assm; ++theta){
			for (i=0; i<=N_T; ++i){	
			rEXT2[theta][i] += difrEXT2(theta,i); 
			}
		}
	} 

	fclose(gaba0_01);
	fclose(gaba0_11);
	fclose(gaba0_21);
	fclose(gaba0_31);
	fclose(gaba0_41);
	fclose(gaba0_51);
	fclose(gaba0_61);
	fclose(gaba0_71);
	fclose(gaba0_81);
	fclose(gaba0_91);
	fclose(gaba1_01);
	fclose(gaba1_11);
	fclose(gaba1_21);
	fclose(gaba1_31);
	fclose(gaba1_41);
	fclose(gaba1_51);
	fclose(gaba1_61);
	fclose(gaba1_71);
	fclose(gaba1_81);
	fclose(gaba1_91);
	fclose(gaba2_01);
	fclose(gaba2_11);
	fclose(gaba2_21);
	fclose(gaba2_31);
	fclose(gaba2_41);
	fclose(gaba2_51);
	fclose(gaba2_61);
	fclose(gaba2_71);
	fclose(gaba2_81);
	fclose(gaba2_91);
	fclose(gaba3_01);
	fclose(gaba3_11);
	fclose(gaba3_21);
	fclose(gaba3_31);
	fclose(gaba3_41);
	fclose(gaba3_51);
	fclose(gaba3_61);
	fclose(gaba3_71);
	fclose(gaba3_81);
	fclose(gaba3_91);
	fclose(gaba4_01);
	fclose(gaba4_11);
	fclose(gaba4_21);
	fclose(gaba4_31);
	fclose(gaba4_41);
	fclose(gaba4_51);
	fclose(gaba4_61);
	fclose(gaba4_71);
	fclose(gaba4_81);
	fclose(gaba4_91);
	fclose(gaba5_01);
	fclose(gaba5_11);
	fclose(gaba5_21);
	fclose(gaba5_31);
	fclose(gaba5_41);
	fclose(gaba5_51);
	fclose(gaba5_61);
	fclose(gaba5_71);
	fclose(gaba5_81);
	fclose(gaba5_91);
	fclose(gaba6_01);
	fclose(gaba6_11);
	fclose(gaba6_21);
	fclose(gaba6_31);
	fclose(gaba6_41);
	fclose(gaba6_51);
	fclose(gaba6_61);
	fclose(gaba6_71);
	fclose(gaba6_81);
	fclose(gaba6_91);
	fclose(gaba7_01);
	fclose(gaba7_11);
	fclose(gaba7_21);
	fclose(gaba7_31);
	fclose(gaba7_41);
	fclose(gaba7_51);
	fclose(gaba7_61);
	fclose(gaba7_71);
	fclose(gaba7_81);
	fclose(gaba7_91);
    
	fclose(uPY0_01);
	fclose(uPY0_11);
	fclose(uPY0_21);
	fclose(uPY0_31);
	fclose(uPY0_41);
	fclose(uPY0_51);
	fclose(uPY0_61);
	fclose(uPY0_71);
	fclose(uPY0_81);
	fclose(uPY0_91);
	fclose(uPY1_01);
	fclose(uPY1_11);
	fclose(uPY1_21);
	fclose(uPY1_31);
	fclose(uPY1_41);
	fclose(uPY1_51);
	fclose(uPY1_61);
	fclose(uPY1_71);
	fclose(uPY1_81);
	fclose(uPY1_91);
	fclose(uPY2_01);
	fclose(uPY2_11);
	fclose(uPY2_21);
	fclose(uPY2_31);
	fclose(uPY2_41);
	fclose(uPY2_51);
	fclose(uPY2_61);
	fclose(uPY2_71);
	fclose(uPY2_81);
	fclose(uPY2_91);
	fclose(uPY3_01);
	fclose(uPY3_11);
	fclose(uPY3_21);
	fclose(uPY3_31);
	fclose(uPY3_41);
	fclose(uPY3_51);
	fclose(uPY3_61);
	fclose(uPY3_71);
	fclose(uPY3_81);
	fclose(uPY3_91);
	fclose(uPY4_01);
	fclose(uPY4_11);
	fclose(uPY4_21);
	fclose(uPY4_31);
	fclose(uPY4_41);
	fclose(uPY4_51);
	fclose(uPY4_61);
	fclose(uPY4_71);
	fclose(uPY4_81);
	fclose(uPY4_91);
	fclose(uPY5_01);
	fclose(uPY5_11);
	fclose(uPY5_21);
	fclose(uPY5_31);
	fclose(uPY5_41);
	fclose(uPY5_51);
	fclose(uPY5_61);
	fclose(uPY5_71);
	fclose(uPY5_81);
	fclose(uPY5_91);
	fclose(uPY6_01);
	fclose(uPY6_11);
	fclose(uPY6_21);
	fclose(uPY6_31);
	fclose(uPY6_41);
	fclose(uPY6_51);
	fclose(uPY6_61);
	fclose(uPY6_71);
	fclose(uPY6_81);
	fclose(uPY6_91);
	fclose(uPY7_01);
	fclose(uPY7_11);
	fclose(uPY7_21);
	fclose(uPY7_31);
	fclose(uPY7_41);
	fclose(uPY7_51);
	fclose(uPY7_61);
	fclose(uPY7_71);
	fclose(uPY7_81);
	fclose(uPY7_91);

	fclose(uLB0_01);
	fclose(uLB1_01);
	fclose(uLB2_01);
	fclose(uLB3_01);
	fclose(uLB4_01);
	fclose(uLB5_01);
	fclose(uLB6_01);
	fclose(uLB7_01);

	fclose(uPY0_02);
	fclose(uPY0_12);
	fclose(uPY0_22);
	fclose(uPY0_32);
	fclose(uPY0_42);
	fclose(uPY0_52);
	fclose(uPY0_62);
	fclose(uPY0_72);
	fclose(uPY0_82);
	fclose(uPY0_92);
	fclose(uPY1_02);
	fclose(uPY1_12);
	fclose(uPY1_22);
	fclose(uPY1_32);
	fclose(uPY1_42);
	fclose(uPY1_52);
	fclose(uPY1_62);
	fclose(uPY1_72);
	fclose(uPY1_82);
	fclose(uPY1_92);
	fclose(uPY2_02);
	fclose(uPY2_12);
	fclose(uPY2_22);
	fclose(uPY2_32);
	fclose(uPY2_42);
	fclose(uPY2_52);
	fclose(uPY2_62);
	fclose(uPY2_72);
	fclose(uPY2_82);
	fclose(uPY2_92);
	fclose(uPY3_02);
	fclose(uPY3_12);
	fclose(uPY3_22);
	fclose(uPY3_32);
	fclose(uPY3_42);
	fclose(uPY3_52);
	fclose(uPY3_62);
	fclose(uPY3_72);
	fclose(uPY3_82);
	fclose(uPY3_92);
	
	fclose(uPY4_02);
	fclose(uPY4_12);
	fclose(uPY4_22);
	fclose(uPY4_32);
	fclose(uPY4_42);
	fclose(uPY4_52);
	fclose(uPY4_62);
	fclose(uPY4_72);
	fclose(uPY4_82);
	fclose(uPY4_92);
	fclose(uPY5_02);
	fclose(uPY5_12);
	fclose(uPY5_22);
	fclose(uPY5_32);
	fclose(uPY5_42);
	fclose(uPY5_52);
	fclose(uPY5_62);
	fclose(uPY5_72);
	fclose(uPY5_82);
	fclose(uPY5_92);
	fclose(uPY6_02);
	fclose(uPY6_12);
	fclose(uPY6_22);
	fclose(uPY6_32);
	fclose(uPY6_42);
	fclose(uPY6_52);
	fclose(uPY6_62);
	fclose(uPY6_72);
	fclose(uPY6_82);
	fclose(uPY6_92);
	fclose(uPY7_02);
	fclose(uPY7_12);
	fclose(uPY7_22);
	fclose(uPY7_32);
	fclose(uPY7_42);
	fclose(uPY7_52);
	fclose(uPY7_62);
	fclose(uPY7_72);
	fclose(uPY7_82);
	fclose(uPY7_92);
	
	fclose(uLB0_02);
	fclose(uLB1_02);
	fclose(uLB2_02);
	fclose(uLB3_02);
	fclose(uLB4_02);
	fclose(uLB5_02);
	fclose(uLB6_02);
	fclose(uLB7_02);

	printf("\a");
    printf("\a");
    printf("\a");
} 

double difuPY1(int thetaa, int ii){
	int jj,thetdash;
	double duPY1=0.0,rec_exT1=0.0,lat_exT=0.0,lat_ihTT=0.0,lat_ihTP=0.0,exP=0.0,c0=0.5,c1=1.0,pw=2.0;
	double top_down = 0.0;
	double alpha_T=1.0; // 1.0 
	double tau_T0=14.0,tau_T1=0.1,tau_T2=0.01,tau_T3=0.01;    

	duPY_leak1[thetaa][ii] = -(DT*gmPY/cmPY)*(uPY1[thetaa][ii]-UPYres); 
	for (jj=0; jj<N_T; ++jj){ 
		rec_exT1 += w_rec_1[thetaa][ii][thetaa][jj]*rPY1[thetaa][jj]; 
	}
	for (jj=0; jj<N_T; ++jj){ 
		rec_exT1 += w_rec_1[thetaa][ii][N_assm][jj]*rPY1[N_assm][jj]; 
	}
	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){
		for (jj=0; jj<=N_T; ++jj){
			if (thetaa!=thetdash) rec_exT1 += w_rec_1[thetaa][ii][thetdash][jj]*rPY1[thetdash][jj];
		}
	}
	duPY_rec_1[thetaa][ii] = -(DT/cmPY)*gAMPA*(uPY1[thetaa][ii]-u_AMPA)*rec_exT1; 
	duPY_fed_1[thetaa][ii] = -(DT/cmPY)*gGABA*(uPY1[thetaa][ii]-u_GABA)*w_fed_1[thetaa][ii]*rSB1[thetaa][ii]; 
	for (jj=0; jj<N_T; ++jj){ 
		lat_ihTT += w_lat_1[thetaa][ii][jj]*rLB1[thetaa][jj];
	}
	duPY_lat_1[thetaa][ii] = -(DT/cmPY)*gGABA*(uPY1[thetaa][ii]-u_GABA)*lat_ihTT; 
	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){
		for (jj=0; jj<=N_T; ++jj){
			top_down += w_v1P_v2[thetaa][ii][thetdash][jj]*rPY2d[thetdash][jj][t-delay];
		}
	}
	duPY_topdown[thetaa][ii] = -(DT/cmPY)*gAMPA*(uPY1[thetaa][ii]-u_AMPA)*top_down; 
	duPY_ext_1[thetaa][ii] = -(DT/cmPY)*gGABA*(uPY1[thetaa][ii]-u_GABA)*delta_T*rEXT[thetaa][ii]; 
	if(t>=onset_0 && t<onset_0+period_0 && thetaa<N_assm){ 
		I_MGB1[thetaa]  = int_inp0_0*exp(-abs(thetaa-theta_inp0)/tau_T0); 
		I_MGB1[thetaa] += int_inp0_1*exp(-abs(thetaa-theta_inp1)/tau_T0); 
		I_MGB1[thetaa] += int_inp0_2*exp(-abs(thetaa-theta_inp2)/tau_T0); 
		I_MGB1[thetaa] += int_inp0_3*exp(-abs(thetaa-theta_inp3)/tau_T0); 
		I_MGB1[thetaa] += int_inp0_4*exp(-abs(thetaa-theta_inp4)/tau_T0); 
		I_MGB1[thetaa] += int_inp0_5*exp(-abs(thetaa-theta_inp5)/tau_T0); 
		I_MGB1[thetaa] += int_inp0_6*exp(-abs(thetaa-theta_inp6)/tau_T0); 
		I_MGB1[thetaa] += int_inp0_7*exp(-abs(thetaa-theta_inp7)/tau_T0); 
	}else{
		I_MGB1[thetaa] = 0.0;
	}
	duPY_MGB1[thetaa] = (DT/cmPY)*(I_MGB1[thetaa]); 
	if (ii>N_T) duPY_MGB1[thetaa] = 0.0; 
	duPY1 = duPY_leak1[thetaa][ii]               
		    + duPY_rec_1[thetaa][ii] 
	        + duPY_fed_1[thetaa][ii] 
		    + duPY_lat_1[thetaa][ii]  
	        + duPY_topdown[thetaa][ii]
		    + duPY_ext_1[thetaa][ii] 
		    + duPY_MGB1[thetaa]*(int)(rand()/32768.0+inp_prob); 
	return(duPY1);
}

double difuPY2(int thetaa, int ii){
	int jj,thetdash;
	double duPY2=0.0,rec_exT2=0.0,lat_exT=0.0,lat_ihTT=0.0,lat_ihTP=0.0,exP=0.0,c0=0.5,c1=1.0,pw=2.0;
	double alpha_T=1.0, tau_T=1.0;   

	duPY_leak2[thetaa][ii] = -(DT*gmPY/cmPY)*(uPY2[thetaa][ii]-UPYres); 
	for (jj=0; jj<N_T; ++jj){ 
		rec_exT2 += w_rec_2[thetaa][ii][thetaa][jj]*rPY2[thetaa][jj]; 
	}
	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){
		for (jj=0; jj<=N_T; ++jj){
			rec_exT2 += w_v2_v1[thetaa][ii][thetdash][jj]*rPY1d[thetdash][jj][t-delay];
		}
	}
	duPY_rec_2[thetaa][ii] = -(DT/cmPY)*gAMPA*(uPY2[thetaa][ii]-u_AMPA)*rec_exT2; 
	for (jj=0; jj<N_T; ++jj){ 
		lat_ihTT += w_lat_2[thetaa][ii][jj]*rLB2[thetaa][jj];
	}
	duPY_lat_2[thetaa][ii] = -(DT/cmPY)*gGABA*(uPY2[thetaa][ii]-u_GABA)*lat_ihTT; 

	duPY_ext_2[thetaa][ii] = -(DT/cmPY)*gGABA*(uPY2[thetaa][ii]-u_GABA)*delta_T2*rEXT2[thetaa][ii]; 
	duPY2 = duPY_leak2[thetaa][ii]                   
		 + duPY_rec_2[thetaa][ii] 
		 + duPY_fed_2[thetaa][ii] 
		 + duPY_lat_2[thetaa][ii]  
		 + duPY_ext_2[thetaa][ii]; 
	return(duPY2);
}

double difuSB1(int thetaa, int ii){
	double duSB1; 
	int thetdash,jj;
	double top_down = 0.0;

	duSB_leak1[thetaa][ii] = -(DT*gmSB/cmSB)*(uSB1[thetaa][ii]-USBres); 
	if (thetaa<N_assm){
		duSB_1[thetaa][ii] = -(DT/cmSB)*gAMPA*(uSB1[thetaa][ii]-u_AMPA)
						   *(wSB_PY1[thetaa][ii]*rPY1[thetaa][ii]); 
	}
	if (thetaa==N_assm){
		duSB_1[thetaa][ii] = -(DT/cmSB)*gAMPA*(uSB1[thetaa][ii]-u_AMPA)
						   *(wSB_PY1[thetaa][ii]*rPY1[thetaa][ii]); 
	}
	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){
		for (jj=0; jj<=N_T; ++jj){
			top_down += w_v1Ia_v2[thetaa][ii][thetdash][jj]*rPY2[thetdash][jj];
		}
	}
	duSB_topdown[thetaa][ii] = -(DT/cmSB)*gAMPA*(uSB1[thetaa][ii]-u_AMPA)*top_down; 
	duSB_ext_1 = -(DT/cmSB)*gGABA*(uSB1[thetaa][ii]-u_GABA)*delta_SB*rEXT[thetaa][ii]; 
	duSB1 = duSB_leak1[thetaa][ii] + duSB_1[thetaa][ii] + duSB_topdown[thetaa][ii] + duSB_ext_1;
	return(duSB1);
}

double difuSB2(int thetaa, int ii){
	double duSB2; 
	duSB_leak2[thetaa][ii] = -(DT*gmSB/cmSB)*(uSB2[thetaa][ii]-USBres);
    
	if (thetaa<N_assm){
		duSB_2[thetaa][ii] = -(DT/cmSB)*gAMPA*(uSB2[thetaa][ii]-u_AMPA)*(wSB_PY1[thetaa][ii]*rPY2[thetaa][ii]); 
	}
	if (thetaa==N_assm){
		duSB_2[thetaa][ii] = -(DT/cmSB)*gAMPA*(uSB2[thetaa][ii]-u_AMPA)*(wSB_PY1[thetaa][ii]*rPY2[thetaa][ii]); 
	}
	duSB_ext_2 = -(DT/cmSB)*gGABA*(uSB2[thetaa][ii]-u_GABA)*delta_SB*rEXT[thetaa][ii]; 
	duSB2 = duSB_leak2[thetaa][ii] + duSB_2[thetaa][ii] + duSB_ext_2;
	return(duSB2);
}

double difuLB1(int thetaa, int ii){
	int jj,thetdash;
	double duLBB,lat_ih=0.0;
	double top_down = 0.0;
	
	duLB_leak1[thetaa][ii] = -(DT*gmLB/cmLB)*(uLB1[thetaa][ii]-ULBres); 
	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){  
		for (jj=0; jj<=N_T; ++jj){
			lat_ih += (wLB_PY1[thetaa][thetdash][ii]*rPY1[thetdash][jj] + wLB_PY1[thetaa][N_assm][ii]*rPY1[N_assm][jj]); 
		}
	}
	duLB1[thetaa][ii] = -(DT/cmLB)*gAMPA*(uLB1[thetaa][ii]-u_AMPA)*lat_ih; 
	duLBB = duLB_leak1[thetaa][ii] + duLB1[thetaa][ii] + duLB1_topdown[thetaa][ii] + duLB_ext1;
	return(duLBB);
}

double difuLB2(int thetaa, int ii){
	int jj,thetdash;
	double duLBB,lat_ih=0.0,bottom_up=0.0;
	
	duLB_leak2[thetaa][ii] = -(DT*gmLB/cmLB)*(uLB2[thetaa][ii]-ULBres); 
	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){  
		for (jj=0; jj<=N_T; ++jj){
			lat_ih += (wLB_PY2[thetaa][thetdash][ii]*rPY2[thetdash][jj] + wLB_PY2[thetaa][N_assm][ii]*rPY2[N_assm][jj]); 
		}
	}
	duLB2[thetaa][ii] = -(DT/cmLB)*gAMPA*(uLB2[thetaa][ii]-u_AMPA)*lat_ih; 
	duLBB = duLB_leak2[thetaa][ii] + duLB2[thetaa][ii] + duLB2_bottom_up[thetaa][ii] + duLB_ext2;
	return(duLBB);
}

double difuG(int thetaa, int ii){ 
	int jj,thetdash;
	double duG=0.0,IatoG=0.0,GtoG=0.0;
	double top_down = 0.0;

	duG_leak[thetaa][ii] = -(DT*gmG/cmG)*(uG[thetaa][ii]-UGres); 
	if (thetaa == 0){
		duG_P[thetaa][ii] = -(DT/cmG)*gAMPA*(uG[thetaa][ii]-u_AMPA)
						   *(wG_P[thetaa]*rPY1[1][ii]
						    +wG_P[thetaa]*rPY1[2][ii]
						    +wG_P[thetaa]*rPY1[3][ii]
						    +wG_P[thetaa]*rPY1[4][ii]
						    +wG_P[thetaa]*rPY1[5][ii]
						    +wG_P[thetaa]*rPY1[6][ii]
						    +wG_P[thetaa]*rPY1[7][ii]
							); 
	}
	if (thetaa == 1){
		duG_P[thetaa][ii] = -(DT/cmG)*gAMPA*(uG[thetaa][ii]-u_AMPA)
						   *(wG_P[thetaa]*rPY1[0][ii]
						    +wG_P[thetaa]*rPY1[2][ii]
						    +wG_P[thetaa]*rPY1[3][ii]
						    +wG_P[thetaa]*rPY1[4][ii]
						    +wG_P[thetaa]*rPY1[5][ii]
						    +wG_P[thetaa]*rPY1[6][ii]
						    +wG_P[thetaa]*rPY1[7][ii]
							); 
	}
	if (thetaa == 2){
		duG_P[thetaa][ii] = -(DT/cmG)*gAMPA*(uG[thetaa][ii]-u_AMPA)
						   *(wG_P[thetaa]*rPY1[0][ii]
						    +wG_P[thetaa]*rPY1[1][ii]
						    +wG_P[thetaa]*rPY1[3][ii]
						    +wG_P[thetaa]*rPY1[4][ii]
						    +wG_P[thetaa]*rPY1[5][ii]
						    +wG_P[thetaa]*rPY1[6][ii]
						    +wG_P[thetaa]*rPY1[7][ii]
							); 
	}
	if (thetaa == 3){
		duG_P[thetaa][ii] = -(DT/cmG)*gAMPA*(uG[thetaa][ii]-u_AMPA)
						   *(wG_P[thetaa]*rPY1[0][ii]
						    +wG_P[thetaa]*rPY1[1][ii]
						    +wG_P[thetaa]*rPY1[2][ii]
						    +wG_P[thetaa]*rPY1[4][ii]
						    +wG_P[thetaa]*rPY1[5][ii]
						    +wG_P[thetaa]*rPY1[6][ii]
						    +wG_P[thetaa]*rPY1[7][ii]
							);
	}
	if (thetaa == 4){
		duG_P[thetaa][ii] = -(DT/cmG)*gAMPA*(uG[thetaa][ii]-u_AMPA)
						   *(wG_P[thetaa]*rPY1[0][ii]
						    +wG_P[thetaa]*rPY1[1][ii]
						    +wG_P[thetaa]*rPY1[2][ii]
						    +wG_P[thetaa]*rPY1[3][ii]
						    +wG_P[thetaa]*rPY1[5][ii]
						    +wG_P[thetaa]*rPY1[6][ii]
						    +wG_P[thetaa]*rPY1[7][ii]
							);
	}
	if (thetaa == 5){
		duG_P[thetaa][ii] = -(DT/cmG)*gAMPA*(uG[thetaa][ii]-u_AMPA)
						   *(wG_P[thetaa]*rPY1[0][ii]
						    +wG_P[thetaa]*rPY1[1][ii]
						    +wG_P[thetaa]*rPY1[2][ii]
						    +wG_P[thetaa]*rPY1[3][ii]
						    +wG_P[thetaa]*rPY1[4][ii]
						    +wG_P[thetaa]*rPY1[6][ii]
						    +wG_P[thetaa]*rPY1[7][ii]
							); 
	}
	if (thetaa == 6){
		duG_P[thetaa][ii] = -(DT/cmG)*gAMPA*(uG[thetaa][ii]-u_AMPA)
						   *(wG_P[thetaa]*rPY1[0][ii]
						    +wG_P[thetaa]*rPY1[1][ii]
						    +wG_P[thetaa]*rPY1[2][ii]
						    +wG_P[thetaa]*rPY1[3][ii]
						    +wG_P[thetaa]*rPY1[4][ii]
						    +wG_P[thetaa]*rPY1[5][ii]
						    +wG_P[thetaa]*rPY1[7][ii]
							); 
	}
	if (thetaa == 7){
		duG_P[thetaa][ii] = -(DT/cmG)*gAMPA*(uG[thetaa][ii]-u_AMPA)
						   *(wG_P[thetaa]*rPY1[0][ii]
						    +wG_P[thetaa]*rPY1[1][ii]
						    +wG_P[thetaa]*rPY1[2][ii]
						    +wG_P[thetaa]*rPY1[3][ii]
						    +wG_P[thetaa]*rPY1[4][ii]
						    +wG_P[thetaa]*rPY1[5][ii]
						    +wG_P[thetaa]*rPY1[6][ii]
							);
	}
	IatoG = wG_Ia[thetaa]*(sF[thetaa][ii]*sF[thetaa][ii]*sF[thetaa][ii]*sF[thetaa][ii])/(sF[thetaa][ii]*sF[thetaa][ii]*sF[thetaa][ii]*sF[thetaa][ii]+Kd); 
	duG_Ia[thetaa][ii] = -(DT/cmG)*gGABAb*(uG[thetaa][ii]-u_GABAb)*IatoG; 
	for (thetdash=0; thetdash<=N_assm-1; ++thetdash){
		for (jj=0; jj<=N_T; ++jj){
			top_down += w_v1G_v2[thetaa][ii][thetdash][jj]*rPY2[thetdash][jj];
		}
	}		
	duG_topdown[thetaa][ii] = -(DT/cmG)*gAMPA*(uG[thetaa][ii]-u_AMPA)*top_down; 
	if (ii == 0){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][0]-uG[thetaa][19]) + (uG[thetaa][0]-uG[thetaa][1])); 
	}
	if (ii == 1){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][1]-uG[thetaa][0]) + (uG[thetaa][1]-uG[thetaa][2])); 
	}
	if (ii == 2){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][2]-uG[thetaa][1]) + (uG[thetaa][2]-uG[thetaa][3])); 
	}
	if (ii == 3){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][3]-uG[thetaa][2]) + (uG[thetaa][3]-uG[thetaa][4])); 
	}
	if (ii == 4){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][4]-uG[thetaa][3]) + (uG[thetaa][4]-uG[thetaa][5])); 
	}
	if (ii == 5){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][5]-uG[thetaa][4]) + (uG[thetaa][5]-uG[thetaa][6])); 
	}
	if (ii == 6){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][6]-uG[thetaa][5]) + (uG[thetaa][6]-uG[thetaa][7])); 
	}
	if (ii == 7){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][7]-uG[thetaa][6]) + (uG[thetaa][7]-uG[thetaa][8])); 
	}
	if (ii == 8){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][8]-uG[thetaa][9]) + (uG[thetaa][8]-uG[thetaa][10])); 
	}
	if (ii == 9){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][9]-uG[thetaa][8]) + (uG[thetaa][9]-uG[thetaa][10])); 
	}
	if (ii == 10){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][10]-uG[thetaa][9]) + (uG[thetaa][10]-uG[thetaa][11])); 
	}
	if (ii == 11){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][11]-uG[thetaa][10]) 
			                + (uG[thetaa][11]-uG[thetaa][12])); 
	}
	if (ii == 12){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][12]-uG[thetaa][11]) + (uG[thetaa][12]-uG[thetaa][13])); 
	}
	if (ii == 13){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][13]-uG[thetaa][12]) + (uG[thetaa][13]-uG[thetaa][14])); 
	}
	if (ii == 14){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][14]-uG[thetaa][13]) + (uG[thetaa][14]-uG[thetaa][15])); 
	}
	if (ii == 15){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][15]-uG[thetaa][14]) + (uG[thetaa][15]-uG[thetaa][16])); 
	}
	if (ii == 16){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][16]-uG[thetaa][15]) + (uG[thetaa][16]-uG[thetaa][17])); 
	}
	if (ii == 17){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][17]-uG[thetaa][16]) + (uG[thetaa][17]-uG[thetaa][18])); 
	}
	if (ii == 18){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][18]-uG[thetaa][17]) + (uG[thetaa][18]-uG[thetaa][19])); 
	}
	if (ii == 19){
		duG_G[thetaa][ii] = -(DT/cmG)*gGap*((uG[thetaa][19]-uG[thetaa][18]) + (uG[thetaa][19]-uG[thetaa][0])); 
	}
	duG = duG_leak[thetaa][ii] + duG_P[thetaa][ii] + duG_Ia[thetaa][ii] + duG_topdown[thetaa][ii] + duG_G[thetaa][ii];
	return(duG);
}

double difrPY1(int thetaa, int ii){
	double drPYT;
	drPYT = DT*(alph_AMPA*GlutPY_c1[thetaa][ii]*(1.0-rPY1[thetaa][ii])-beta_AMPA*rPY1[thetaa][ii]); 
	return(drPYT);
}

double difrPY2(int thetaa, int ii){
	double drPYT;
	drPYT = DT*(alph_AMPA*GlutPY_c2[thetaa][ii]*(1.0-rPY2[thetaa][ii])-beta_AMPA*rPY2[thetaa][ii]); 
	return(drPYT);
}

double difrSB1(int thetaa, int ii){
	double drSBT;
	drSBT = DT*(K1*GABASB_c1[thetaa][ii][t]*(1.0-rSB1[thetaa][ii])-K2*rSB1[thetaa][ii]); 
	return(drSBT);
}

double difsF(int thetaa, int ii){
	double dsF;
	dsF = DT*(K3*rSB1[thetaa][ii]-K4*sF[thetaa][ii]); 
	return(dsF);
}

double difrSB2(int thetaa, int ii){
	double drSBT;
	drSBT = DT*(alph_GABA*GABASB_c2[thetaa][ii][t]*(1.0-rSB2[thetaa][ii])-beta_GABA*rSB2[thetaa][ii]); 
	return(drSBT);
}

double difrLB1(int thetaa, int ii){
	double drLB;
	drLB = DT*(alph_GABA*GABALB_c1[thetaa][ii][t]*(1.0-rLB1[thetaa][ii])-beta_GABA*rLB1[thetaa][ii]); 
	return(drLB);
}

double difrLB2(int thetaa, int ii){
	double drLB;
	drLB = DT*(alph_GABA*GABALB_c2[thetaa][ii][t]*(1.0-rLB2[thetaa][ii])-beta_GABA*rLB2[thetaa][ii]); 
	return(drLB);
}

double difrEXT(int thetaa, int ii){
	double drEXT;
	drEXT = DT*(alph_GABA*GABA_ext[thetaa][ii]*(1.0-rEXT[thetaa][ii])-beta_GABA*rEXT[thetaa][ii]); 
	return(drEXT);
}

double difrEXT2(int thetaa, int ii){
	double drEXT2;
	drEXT2 = DT*(alph_GABA*GABA_V2*(1.0-rEXT2[thetaa][ii])-beta_GABA*rEXT2[thetaa][ii]); 
	return(drEXT2);
}

void dfsGABA_ext(void){
	int theta,i;

	for (theta=0; theta<=N_assm; ++theta){
			for (i=0; i<=N_T; ++i){	
				I_GABA[theta][i] = 0.0;
			}
	}
	for (theta=0; theta<=N_assm; ++theta){
			for (i=0; i<=N_T; ++i){	
				I_GABA[theta][i] += -gamma*(GABA_ext[theta][i]-GABA_c0)*DT 
					                + m_G*(uG[theta][i]-uG_trn)*(GABAamb_max-GABA_ext[theta][i])*(GABA_ext[theta][i]-GABAamb_min)*DT; 
				GABA_extw[theta][i] = GABA_ext[theta][i] + I_GABA[theta][i];
			}
	}
	for (theta=0; theta<=N_assm; ++theta){
			for (i=0; i<=N_T; ++i){	
				GABA_ext[theta][i] = GABA_extw[theta][i];
			}
	}
}

double sigmoidPY(double uPYY){
	return(1/(1+exp(-steep_PY*(uPYY-thres_PY))));
}

double sigmoidPY2(double uPYY){
	return(1/(1+exp(-steep_PY2*(uPYY-thres_PY2))));
}

double sigmoidSB(double uSBB){
	return(1/(1+exp(-steep_SB*(uSBB-thres_SB))));
}

double sigmoidSB2(double uSBB){
	return(1/(1+exp(-steep_SB2*(uSBB-thres_SB2))));
}

double sigmoidLB(double uLBB){
	return(1/(1+exp(-steep_LB*(uLBB-thres_LB))));
}

double sigmoidLB2(double uLBB){
	return(1/(1+exp(-steep_LB2*(uLBB-thres_LB2))));
}

double rand01(long int *ix){
   double x;
   long int of;

   *ix=(*ix)*48828125;
   if(*ix<0){
	   of=2147483647;
       *ix=(*ix+of)+1;
   }
   x=(double)(*ix)*0.4656613e-9;
   return(x);
 }

void init(void){
	int theta,i,thetaa,ii,thetdash,jj;
	double tau_lat=1000.0;   // distribution of P-to-B projectionk
	double w_lat_excit=0.0;  // non-use 
    double alph_w_L_TV1=1.2; // P-to-B connection weight in N_S    
    double alph_w_L_TV2=1.2; // in N_M 
	double w_v2_1=4.6;       // P(N_S)-to-P(N_M) connection weight
	double w_v1Ia_2=0.0;     // non-use
	double w_v1Ib_2=0.0;     // non-use
	double w_v1P_2=4.6;      // 4.6 V2(P) to V1(P)  connection strength
	double w_v1G_2=0.0;      // non-use

	srand(17);
	for (thetaa=0; thetaa<=N_assm; ++thetaa){
			for (ii=0; ii<=N_T; ++ii){
     			for (thetdash=0; thetdash<=N_assm; ++thetdash){
					for (jj=0; jj<=N_T; ++jj){
						w_rec_1[N_assm][ii][thetdash][jj] = 0.0;  
						w_rec_2[N_assm][ii][thetdash][jj] = 0.0;  
					}
				}
			}
	}
	for (thetaa=0; thetaa<=N_assm; ++thetaa){
			for (ii=0; ii<=N_T; ++ii){
     			for (thetdash=0; thetdash<=N_assm; ++thetdash){
					for (jj=0; jj<=N_T; ++jj){
						w_v2_v1[thetaa][ii][thetdash][jj] = 0.0;  
					}
				}
			}
	}
	for (thetaa=0; thetaa<=N_assm; ++thetaa){
			for (ii=0; ii<=N_T; ++ii){
     			for (thetdash=0; thetdash<=N_assm; ++thetdash){
					for (jj=0; jj<=N_T; ++jj){
						w_v1Ia_v2[thetaa][ii][thetdash][jj] = 0.0;  
					}
				}
			}
	}
	for (thetaa=0; thetaa<=N_assm; ++thetaa){
			for (ii=0; ii<=N_T; ++ii){
     			for (thetdash=0; thetdash<=N_assm; ++thetdash){
					for (jj=0; jj<=N_T; ++jj){
						w_v1G_v2[thetaa][ii][thetdash][jj] = 0.0; 
					}
				}
			}
	}
	for (thetaa=4; thetaa<=N_assm-1; ++thetaa){
		for (ii=0; ii<=N_T; ++ii){
     			for (thetdash=0; thetdash<=3; ++thetdash){
					for (jj=0; jj<=N_T; ++jj){
						if (abs(thetaa-thetdash) == 4) 
						{
							w_rec_1[thetaa][ii][thetdash][jj] = w_lat_excit;
						}
					}
				}
		}
	}
	for (thetaa=0; thetaa<=N_assm-1; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			for (jj=0; jj<=N_T; ++jj){
				if (ii!=jj)	w_rec_1[thetaa][ii][thetaa][jj]  = wLPPV1;   
				if (ii!=jj)	w_rec_2[thetaa][ii][thetaa][jj]  = wLPPV2; 
			}
		}
	}
	for (thetaa=0; thetaa<=N_assm; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			w_fed_1[thetaa][ii] = 0.0;  
			w_fed_2[thetaa][ii] = 0.0; 
		}
	}
	for (thetaa=0; thetaa<=N_assm-1; ++thetaa){   
		for (ii=0; ii<=N_T; ++ii){
			for (jj=0; jj<=N_T; ++jj){
				w_lat_1[thetaa][ii][jj] = 6.0;    
				w_lat_2[thetaa][ii][jj] = 6.0;   
			}
		}
	}
	for (thetaa=0; thetaa<=N_assm-1; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			wSB_PY1[thetaa][ii] = 60.0; 
			wSB_PY2[thetaa][ii] = 0.0;   
		}
	}
	for (ii=0; ii<=N_T; ++ii){
		wLB_PY1[0][5][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[0][6][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0)));  
		wLB_PY1[0][7][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[0][1][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[0][2][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[0][3][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[0][4][ii] = alph_w_L_TV1*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY1[1][6][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[1][7][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[1][0][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[1][2][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[1][3][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[1][4][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0)));  
		wLB_PY1[1][5][ii] = alph_w_L_TV1*exp(-(pow((0-4)/tau_lat,2.0)));   

		wLB_PY1[2][7][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[2][0][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[2][1][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[2][3][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[2][4][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[2][5][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0)));  
		wLB_PY1[2][6][ii] = alph_w_L_TV1*exp(-(pow((0-4)/tau_lat,2.0)));  

		wLB_PY1[3][0][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[3][1][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[3][2][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[3][4][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0)));  
		wLB_PY1[3][5][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0)));  
		wLB_PY1[3][6][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[3][7][ii] = alph_w_L_TV1*exp(-(pow((0-4)/tau_lat,2.0)));  

		wLB_PY1[4][1][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[4][2][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[4][3][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[4][5][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[4][6][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[4][7][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[4][0][ii] = alph_w_L_TV1*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY1[5][2][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[5][3][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0)));  
		wLB_PY1[5][4][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[5][6][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[5][7][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[5][0][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[5][1][ii] = alph_w_L_TV1*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY1[6][3][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[6][4][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[6][5][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[6][7][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[6][0][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[6][1][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[6][2][ii] = alph_w_L_TV1*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY1[7][4][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY1[7][5][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY1[7][6][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[7][0][ii] = alph_w_L_TV1*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY1[7][1][ii] = alph_w_L_TV1*exp(-(pow((0-2)/tau_lat,2.0)));  
		wLB_PY1[7][2][ii] = alph_w_L_TV1*exp(-(pow((0-3)/tau_lat,2.0)));  
		wLB_PY1[7][3][ii] = alph_w_L_TV1*exp(-(pow((0-4)/tau_lat,2.0)));  

		wLB_PY2[0][5][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[0][6][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[0][7][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[0][1][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[0][2][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[0][3][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[0][4][ii] = alph_w_L_TV2*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY2[1][6][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[1][7][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[1][0][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[1][2][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[1][3][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[1][4][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[1][5][ii] = alph_w_L_TV2*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY2[2][7][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[2][0][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[2][1][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[2][3][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[2][4][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[2][5][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[2][6][ii] = alph_w_L_TV2*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY2[3][0][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[3][1][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[3][2][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[3][4][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[3][5][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[3][6][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[3][7][ii] = alph_w_L_TV2*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY2[4][1][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[4][2][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[4][3][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[4][5][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[4][6][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[4][7][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[4][0][ii] = alph_w_L_TV2*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY2[5][2][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[5][3][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[5][4][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[5][6][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[5][7][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[5][0][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[5][1][ii] = alph_w_L_TV2*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY2[6][3][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[6][4][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[6][5][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[6][7][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[6][0][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[6][1][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[6][2][ii] = alph_w_L_TV2*exp(-(pow((0-4)/tau_lat,2.0))); 

		wLB_PY2[7][4][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[7][5][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[7][6][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[7][0][ii] = alph_w_L_TV2*exp(-(pow((0-1)/tau_lat,2.0))); 
		wLB_PY2[7][1][ii] = alph_w_L_TV2*exp(-(pow((0-2)/tau_lat,2.0))); 
		wLB_PY2[7][2][ii] = alph_w_L_TV2*exp(-(pow((0-3)/tau_lat,2.0))); 
		wLB_PY2[7][3][ii] = alph_w_L_TV2*exp(-(pow((0-4)/tau_lat,2.0))); 
	}
	for (thetaa=0; thetaa<=N_assm; ++thetaa){  
		wG_P[thetaa]  = 0.0; 
		wG_Ia[thetaa] = 3.0; 
	}
	for (thetaa=0; thetaa<=N_assm-1; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			for (jj=0; jj<=N_T; ++jj) w_v1P_v2[thetaa][ii][thetaa][jj] = w_v1P_2;
		}
	}
	for (thetaa=0; thetaa<=N_assm-1; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			for (jj=0; jj<=N_T; ++jj) w_v2_v1[thetaa][ii][thetaa][jj] = w_v2_1;
		}
	}
	for (theta=0; theta<=N_assm; ++theta){
		for (i=0; i<=N_T; ++i){	
			GABA_ext[theta][i]  = GABA_c0; 
			GABA_extw[theta][i] = GABA_c0; 
		}
	}
}	

void display(void){
	int thetaa,ii,jj;
	
	for (thetaa=0; thetaa<=N_assm; ++thetaa){  
		for (ii=0; ii<=N_T; ++ii){
			for (jj=0; jj<=N_T; ++jj){
				if (ii!=jj){
					if ((thetaa==NEURONw) & (t%INTV==0)){ 
						printf("t=%d w_rec_1[%d][%d][%d][%d]=%lf \n",t,thetaa,ii,thetaa, jj,
						w_rec_1[thetaa][ii][thetaa][jj]);
					}
				}    
			}
		}	
	}
}