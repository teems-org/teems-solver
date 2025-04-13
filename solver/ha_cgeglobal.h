#ifndef HA_CGEGLOBAL_H_INCLUDED
#define HA_CGEGLOBAL_H_INCLUDED

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <petscksp.h>
#include <petscsys.h>
#include <petsctime.h>
#include <stdbool.h>
#include </home/teems-solver/lib/petsc/src/mat/impls/aij/seq/aij.h>
#include <petscdmda.h>
#include <omp.h>

#define NAMESIZE 256
#define DATREADLINE 150000
#define TABREADLINE 20000//2536
#define HEADERSIZE 5
#define TABLINESIZE 20000//81
#define NOPERTINSUM 5
#define MAXVARDIM 10 //maximum variable dimension
#define MAXSUPSET 12 //original
#define SORD 1
#define MAXSSIZE 187500000//1500000000/8
//To switch between single and double precision: SORD 0 single; SORD 1 double; change ha_cgetype and ha_floattype to float or double; change hsl_mp48ss.f90 accordingly. Note should change GUI too!!!
int isLinux;
int smallthreads;
int medthreads;
int mymaxnumthrd;
double kindx1,kindx2,kval1,kval2,kval3;
int step1,step2,step3;
MPI_Comm HA_COMM,HA1_COMM;
static char temdir[]="/tmp/";

typedef double ha_cgetype;
typedef int uvdim;
typedef long int uvadd;
typedef int uvaddshort;//relax if nvarele greater than 2 billions
typedef long int forint;
typedef float ha_floattype;
//typedef int utlin;
//typedef int svdim;
//typedef long int svadd;
//typedef int stlin;
//typedef long int uvar;

//thu
//typedef struct
//{
//    char varname[NAMESIZE];
//    char varindex[NAMESIZE];
//    char varset[NAMESIZE];
//    char varpos[NAMESIZE];
//    ha_cgetype varval;
//} ha_cgevarname;

//ha_cgeiof.c
typedef struct
{
  char logname[NAMESIZE];
  char filname[TABREADLINE];
} hcge_iodata ;
int hcge_niodata(char *fname,char *comsyntax);

typedef struct
{
  uvdim dim1;
  char ch[NAMESIZE];
} ha_cgemvar1 ;
int ha_cgermvar1(char *varname, char *filename,uvdim d1, ha_cgemvar1 *record);
int ha_cgerdvar1(char *varname, char *filename,uvdim *vsize, char *longname,uvdim *d1);
//int ha_cgewmvar1(char *vname, char *longname, char *fname, char *wmode, int rsize, int vsize, int d1, ha_cgemvar1 *record);
int hcge_rcmd(char *filename, int niodata, hcge_iodata *iodata, char *tabfile, char *closure, char *shock);
int hcge_wtab(char *filename, char *newtabfile);
int hcge_rsetname(char *filename, char *varname, int indx, char *setname);

//ha_cgetab.c

typedef struct
{
  //char varname[NAMESIZE];
  ha_floattype varval;
} ha_cgesumele ;
typedef struct
{
  //char varname[NAMESIZE];
  ha_floattype varval;
  ha_floattype var0;
  //ha_cgetype varchange;
  ha_floattype csolpupd;//cursol;//In stochastic case, with plus var, var0 holds x value.
  //ha_floattype pastupd;
} ha_cgevar ;
char* ha_revstrpbrk(char *line, char *finditems);
int ha_cgerevfind(char *line, char *finditem);
int ha_cgenchf(char *line, int finditem);
//int ha_cgefpmul(char *fomulain, ha_cgevar *record, int sofrecord, ha_cgevar *varpar, int npar, ha_cgevar *varpow, int npow, ha_cgevar *varmul, int nmul);
int ha_cgenfind(char *line, char *finditem);
char* ha_cgeeqfind(char *line, int finditem);
bool intreadCSV(char *fileName,int* vec, int vecCol);
//int ha_cgefpdiv(char *fomulain, ha_cgevar *record, int sofrecord, ha_cgevar *varpar, int npar, ha_cgevar *varpow, int npow, ha_cgevar *varmul, int nmul, ha_cgevar *vardiv, int ndiv);
//int ha_cgefpplu(char *fomulain, ha_cgevar *record, int sofrecord, ha_cgevar *varpar, int npar, ha_cgevar *varpow, int npow, ha_cgevar *varmul, int nmul, ha_cgevar *varplu, int nplu);
//int ha_cgefpmin(char *fomulain, ha_cgevar *record, int sofrecord, ha_cgevar *varpar, int npar, ha_cgevar *varpow, int npow, ha_cgevar *varmul, int nmul, ha_cgevar *varplu, int nplu, ha_cgevar *varmin, int nmin);

//ha_cgetab.c

typedef struct
{
  char header[HEADERSIZE];
  int fileid;
  char setname[NAMESIZE];
  char readele[TABREADLINE];
  uvadd begadd;
  uvdim size;
  uvdim subsetid[MAXSUPSET];//Supersetid to be more precised
  //char supersetname[MAXSUPSET][NAMESIZE];
  //uvdim supersetsize[MAXSUPSET];
  bool intertemp;
  int intsup;
  bool regional;
  int regsup;
} ha_cgeset ;
typedef struct
{
  char setele[NAMESIZE];
  uvdim setsh[MAXSUPSET];
} ha_cgesetele ;

typedef struct
{
  char cofname[NAMESIZE];
  uvadd begadd;
  uvdim size;
  //uvdim dims[MAXVARDIM];//dimension size
  //uvadd setbegadd;
  uvadd setid[MAXVARDIM];
  //uvadd setbegadd1[MAXVARDIM];
  uvadd antidims[MAXVARDIM];
  //char dimsets[MAXVARDIM][NAMESIZE];
  uvadd matsize;
  bool level_par;
  bool change_real;
  bool suplval;
  int gltype;//1 GE 2 GT 3 LE 4 LT
  ha_floattype glval;
} hcge_cof ;

typedef struct
{
  char sumname[NAMESIZE];
  //uvdim sumsize;
  uvadd summatsize;
  //char sumset[NAMESIZE];
  char sumindx[NAMESIZE];
  //uvdim sumsubsetid;
  //uvdim sumsupsetsize;
  uvdim size;
  uvadd begadd;
  //uvadd ssetbegadd;
  uvadd sumsetid;
  char dimnames[MAXVARDIM][NAMESIZE];
  uvadd setid[MAXVARDIM];
  //char dimsets[MAXVARDIM][NAMESIZE];
  //uvdim dims[MAXVARDIM];
  //uvadd dimssetbegadd[MAXVARDIM];
  //uvdim subsetid[MAXVARDIM][MAXSUPSET];
  //uvdim supsetsize[MAXVARDIM];
  uvadd antidims[MAXVARDIM];
} hcge_sumcof ;

typedef struct
{
  //char cofname[NAMESIZE];
  ha_floattype cofval;
} ha_cgecofele ;

//typedef struct
//{
//  long int row;
//  long int col;
//  ha_cgetype mval;
//} ha_cgespamat ;

typedef struct
{
  //char ExoName[NAMESIZE];
  uvaddshort ExoIndx;
  bool ShockId;//1 exo
  ha_floattype ShockVal;
} ha_cgeexovar ;

typedef struct
{
  //char arSet[NAMESIZE];
  char arIndx[NAMESIZE];
  uvadd setid;
  uvdim indx;
  //uvdim SetSize;
  //uvadd SetBegAdd;
  //uvdim SuperSetSize;
  //uvdim subsetid[MAXSUPSET];
} ha_cgesetindx ;

//typedef struct
//{
//  char SetIndxName[NAMESIZE];
//  uvadd SetBegAdd;
//  uvdim SetSize;
//  uvdim SetIndxNum;
//} ha_cgesumindx ;

typedef struct
{
  char LinVarName[NAMESIZE];
  uvadd LinVarIndx;
  char dimnames[MAXVARDIM][NAMESIZE];
  char dimsetnames[MAXVARDIM][NAMESIZE];
  uvdim dimleadlag[MAXVARDIM];
  uvdim dimindx[MAXVARDIM];
} hcge_linvars ;

//typedef struct
//{
//  uvadd begadd;
//  uvdim vartype;
//  uvdim size;
//  uvadd antidims[MAXVARDIM];
//  uvdim arSet[MAXVARDIM];
//} hcge_cvar ;
int ha_csumindx(char *formulain);
int hcge_wvar(char *filename, char *newtabfile,hcge_cof *ha_var,uvadd nvar);

uvdim ha_cgenset(char *fname);
int ha_cgerset(char *fname, int niodata, hcge_iodata *iodata, ha_cgeset *record,uvdim nset);
int hcge_rinterset(char *fname, int niodata, hcge_iodata *iodata, ha_cgeset *record,uvdim nset);
uvdim ha_setunion(ha_cgesetele *ha_setele, ha_cgeset *ha_set,uvdim nset,uvdim i);
uvdim ha_setplus(ha_cgesetele *ha_setele, ha_cgeset *ha_set,uvdim nset,uvdim i);
uvdim ha_setminus(ha_cgesetele *ha_setele, ha_cgeset *ha_set,uvdim nset,uvdim i);
uvdim ha_cgesubsetchck(ha_cgesetele *ha_setele, ha_cgeset *ha_set,uvdim nset,uvadd* contin);
char *ha_cgercls(char *commsyntax, FILE *filehandle, char *readline);
char *ha_cgertabl(char *commsyntax, FILE *filehandle, char *readline,uvadd rlinesize);
char *ha_cgertabl1(char *commsyntax, FILE *filehandle, char *readline, ha_cgevar *record, hcge_cof *ha_cof,uvadd ncof,ha_cgetype *zerodivide,uvadd rlinesize);
int ha_cgefind(char *line, char *finditem);
char *ha_cgefrstr(char *line, char *finditem, char *replitem);
int ha_chrfrall(char *line, int finditem, int replitem);
char *ha_cgefrstrvbz(char *line, char *finditem, char *replitem,uvdim nbuffer);
char *ha_cgefrstrvbz1(char *line, char *finditem, char *replitem,uvdim nbuffer);
char *ha_cgefrstr1(char *line, char *finditem, char *replitem);
char *ha_cgedrcmt(char *line, char *token);
int hcge_defvar(char *fname, hcge_cof *record, uvadd ncof);
uvadd hcge_rvar(char *fname, char *commsyntax, hcge_cof *record, uvadd ncof, ha_cgeset *ha_set,uvdim nset);
uvadd ha_cgeralltime(ha_cgeset *ha_set,uvdim nset);
uvadd ha_cgencof(char *fname, char *commsyntax);
uvadd hcge_rexo(char *fname, char *commsyntax,ha_cgeexovar *ha_cgeshock, hcge_cof *ha_var,uvadd nvar,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele);
uvadd hcge_rshock(char *fname, char *commsyntax,ha_cgeexovar *ha_cgeshock,uvadd nvarele, hcge_cof *ha_var,uvadd nvar,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele,uvdim subints);
uvadd hcge_rcof(char *fname, char *commsyntax, hcge_cof *record, uvadd ncof, ha_cgeset *ha_set,uvadd nset);
int hcge_rcofele(hcge_cof *ha_cof,uvadd ncof, ha_cgeset *ha_set,uvdim nset, ha_cgecofele *ha_cofele);
//int hcge_rsumele(hcge_sumcof *ha_cof,uvadd ncof, ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele,ha_cgesumele *ha_cofele);
//int hcge_mrsumele(hcge_sumcof *ha_cof, long int ncof, ha_cgeset *ha_set, int nset, ha_cgesetele *ha_setele,ha_cgevar *ha_cofele);
uvadd hcge_readff(char *fname, int niodata, hcge_iodata *iodata, char *commsyntax,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele,hcge_cof *ha_cof,uvadd ncof, ha_cgecofele *ha_cofele,uvadd ncofele,hcge_cof *ha_var,uvadd nvar, ha_cgecofele *ha_varele,uvadd nvarele);
uvadd rdatconv(hcge_cof *ha_var,uvadd varindx,uvadd loops,ha_cgeset *ha_set);
int hcge_repllin(char *formulain,int linindx);
int hcge_rlinzero(char *formulain,int linindx);
int hcge_nsum(char *formulain, char *commsyntax);
//int hcge_rsum(char *formulain, char *commsyntax, hcge_sumcof *sum_cof);
uvadd hcge_dsum(char *formulain, char *commsyntax, hcge_sumcof *sum_cof,ha_cgesetindx *arSet,ha_cgeset *ha_set,uvdim nset,uvdim fdim,int j);
int ha_cgecutsum(char *formula);
int ha_cgerecovar(char *fomulain);
uvadd ha_cgersubset(char *fname, ha_cgesetele *ha_setele, ha_cgeset *ha_set,uvdim nset);
char *ha_cgefrchr(char *line, int finditem, int replitem);
//int ha_spadd_size(int *irnA, int *jcnA, int *irnB, int *jcnB,int nrow,int *irnC);
//int ha_spadd_oper(int *irnA, int *jcnA, ha_cgetype *valsA, int *irnB, int *jcnB, ha_cgetype *valsB,int nrow,int *irnC, int *jcnC,ha_cgetype *valsC);
//thu
//Moi

typedef struct
{
  uvdim Oper;//0 none; 1 *;2 /;3 +; 4 -;5 ^;71 if =;72 if >; 73 if <; 74 if <>; 75 if <=; 76 if >=;
  //int TmpVarType;//0 temvar;1 par; 2 number
  char TmpVarName[NAMESIZE];
  ha_floattype TmpVarVal;
  uvdim Var1Type;//0 var cof;1 lin; 2 sum;3 lin;4 calvars; 5 number; 6 varchange; 41 id01, 42 abs, 43 loge
  uvadd Var1BegAdd;
  int Var1leadlag[MAXVARDIM];
  //uvadd Var1Size;
  uvadd Var1SupSet[MAXVARDIM];
  uvadd Var1SSIndx[MAXVARDIM];
  //uvadd Var1SetAdd[MAXVARDIM];
  uvadd Var1ADims[MAXVARDIM];
  ha_floattype Var1Val;
  uvdim Var2Type;
  uvadd Var2BegAdd;
  int Var2leadlag[MAXVARDIM];
  //uvadd Var2Size;
  uvadd Var2SupSet[MAXVARDIM];
  uvadd Var2SSIndx[MAXVARDIM];
  //uvadd Var2SetAdd[MAXVARDIM];
  uvadd Var2ADims[MAXVARDIM];
  ha_floattype Var2Val;

  uvdim Var3Type;
  uvadd Var3BegAdd;
  int Var3leadlag[MAXVARDIM];
  //uvadd Var2Size;
  uvadd Var3SupSet[MAXVARDIM];
  uvadd Var3SSIndx[MAXVARDIM];
  //uvadd Var2SetAdd[MAXVARDIM];
  uvadd Var3ADims[MAXVARDIM];
  ha_floattype Var3Val;
} hcge_calvars ;

ha_cgetype hnew_simplrpl(char *var2, ha_cgevar *record, hcge_cof *ha_cof,uvadd ncof);
int hnew_varrepl(char *var2, ha_cgeset *ha_set,hcge_cof *ha_cof,uvadd ncof, hcge_cof *ha_var,uvadd nvar,uvadd ncofele,hcge_sumcof *sum_cof,int totalsum,hcge_calvars *ha_calvar,int ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim,int varindex);
int hnew_intrpl(char *line);
int hnew_arset(char *p,int *leadlag);
int ha_newfparse(char *fomulain, ha_cgeset *ha_set,hcge_cof *ha_cof, uvadd ncof, hcge_cof *ha_var,uvadd nvar,uvadd ncofele,hcge_sumcof *sum_cof,uvdim totalsum,hcge_calvars *ha_calvar,uvdim *ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim);
ha_cgetype ha_newfpcal(ha_cgevar *record, ha_cgeset *ha_set,ha_cgesetele *ha_setele,ha_cgesumele *ha_sumele,hcge_calvars *ha_calvar,int ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim, ha_cgetype zerodivide);
ha_cgetype ha_newfpcal01(ha_cgevar *record, ha_cgeset *ha_set,ha_cgesetele *ha_setele,ha_cgesumele *ha_sumele,hcge_calvars *ha_calvar,int ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim, ha_cgetype zerodivide);
int ha_newfppow(char *fomulain, ha_cgeset *ha_set,int npow,int ipar,hcge_cof *ha_cof,uvadd ncof, hcge_cof *ha_var,uvadd nvar,uvadd ncofele,hcge_sumcof *sum_cof,int totalsum,hcge_calvars *ha_calvar,int *ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim);
int ha_newfpmuldiv(char *fomulain, ha_cgeset *ha_set,int nmul,int ipar,hcge_cof *ha_cof,uvadd ncof, hcge_cof *ha_var,uvadd nvar,uvadd ncofele,hcge_sumcof *sum_cof,int totalsum,hcge_calvars *ha_calvar,int *ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim);
int ha_newfpplumin(char *fomulain, ha_cgeset *ha_set,int nplu,int ipar,hcge_cof *ha_cof,uvadd ncof, hcge_cof *ha_var,uvadd nvar,uvadd ncofele,hcge_sumcof *sum_cof,int totalsum,hcge_calvars *ha_calvar,int *ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim);
int ha_newfpplumin01(char *fomulain, ha_cgeset *ha_set,int nplu,int ipar,hcge_cof *ha_cof,uvadd ncof, hcge_cof *ha_var,uvadd nvar,uvadd ncofele,hcge_sumcof *sum_cof,int totalsum,hcge_calvars *ha_calvar,int *ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim);
int ha_newfpif(char *fomulain, ha_cgeset *ha_set,int nif,int ipar,hcge_cof *ha_cof,uvadd ncof, hcge_cof *ha_var,uvadd nvar,uvadd ncofele,hcge_sumcof *sum_cof,int totalsum,hcge_calvars *ha_calvar,int *ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim);
uvadd hnew_calcff(char *fname, char *commsyntax,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele,bool IsIni);
int hnew_calsum(char *formulain, char *commsyntax,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele,ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele, hcge_cof *ha_cof,uvadd ncof, hcge_cof *ha_var,uvadd nvar,hcge_sumcof *sum_cof,int totalsum,ha_cgesumele *ha_sumele,uvadd nsumele,hcge_calvars *ha_calvar,ha_cgesetindx *arSet1,uvdim fdim,int *sumindx,int j, ha_cgetype zerodivide);
uvadd hnew_update(char *fname,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele);
uvadd hnew_mupdate(char *fname,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele);
uvadd hnew_graggupd(char *fname,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele);
uvadd hnew_gupd(char *fname,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele);
//lin

int hlin_dsum(char *formulain, char *commsyntax, hcge_sumcof *sum_cof,ha_cgesetindx *arSet,ha_cgeset *ha_set,uvdim nset,uvdim fdim,int j);
int hlin_calsum(char *formulain, char *commsyntax,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele,ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar,hcge_sumcof *sum_cof,int totalsum,ha_cgesumele *ha_sumele,uvadd nsumele,hcge_calvars *ha_calvar,ha_cgesetindx *arSet1,uvdim fdim,int *sumindx,int j, ha_cgetype zerodivide);
int hlin_replsum(char *formulain, char *commsyntax,int LinIndx, hcge_linvars *LinVars,hcge_cof *ha_var);
int HaNewMatVal(char *fname, char *commsyntax,ha_cgeset *ha_set,uvadd nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele,ha_cgeexovar *ha_cgeshock,uvadd ndblock,uvadd alltimeset,uvadd allregset,uvadd *ha_eqadd,uvadd *counteq,uvadd nintraeq,Mat A,Mat B);
int hlin_rlinone(char *formulain,hcge_linvars *LinVars,int linindx,hcge_cof *ha_var);
int NewMatreadele(char *fname, char *commsyntax,ha_cgeset *ha_set,uvdim nset,ha_cgesetele *ha_setele,hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar,ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele, uvadd nexo,ha_cgeexovar *ha_cgeshock,uvadd ndblock,uvadd alltimeset,uvadd allregset,bool *ha_eqint,uvadd *ha_eqadd,uvdim *ha_eqtime,uvdim *ha_eqreg,uvadd *counteq,uvadd nintraeq,bool *sbbd_overrid,PetscInt Istart,PetscInt Iend,PetscInt *dnz,PetscInt *dnnz,PetscInt *onz,PetscInt *onnz,PetscInt *dnzB,PetscInt *dnnzB,PetscInt *onzB,PetscInt *onnzB,int nesteddbbd);
int NewMatvarRead(char *fname, char *commsyntax,ha_cgeset *ha_set,uvdim nset,ha_cgesetele *ha_setele,hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar,ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele,ha_cgeexovar *ha_cgeshock,bool *var_inter,hcge_cof *ha_eq,bool *ha_eqint,uvdim *eq_orderintra,uvdim *eq_orderreg,uvadd allregset,uvadd alltimeset,uvdim *orderintra,uvdim *orderreg);
int NestedMatvarRead(char *fname, char *commsyntax,ha_cgeset *ha_set,uvdim nset,ha_cgesetele *ha_setele,hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar,ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele,ha_cgeexovar *ha_cgeshock,bool *var_inter,hcge_cof *ha_eq,bool *ha_eqint,uvdim *eq_orderintra,uvdim *eq_orderreg,uvadd allregset,uvadd alltimeset,uvdim *orderintra,uvdim *orderreg);
bool NDBBD_sol(PetscInt rank, int begmat,int nreg,int * insize,int insizes, Mat **submatCij,Mat **submatBij,ha_cgetype *b,ha_cgetype *sol,bool ifremove,char** fn01,char** fn02, char** fn03);
bool NDBBD_sol_nread(PetscInt rank, int begmat,int nreg,int * insize,int insizes, Mat **submatCij,Mat **submatBij,ha_cgetype *b,ha_cgetype *sol,int** irnereg,int** keepreg,ha_cgetype** valereg,ha_cgetype *cntl,ha_cgetype *rinfo,ha_cgetype *error1,int *icntl,int *info,ha_cgetype *w,int *iw);
bool NDBBD_sol_nread1(PetscInt rank, int begmat,int nreg,int * insize,int insizes, Mat **submatCij,Mat **submatBij,ha_cgetype *b,ha_cgetype *sol,int** irnereg,int** keepreg,ha_cgetype** valereg,ha_cgetype *cntl,ha_cgetype *rinfo,ha_cgetype *error1,int *icntl,int *info,ha_cgetype *w,int *iw,ha_cgetype *b02);
uvadd hcge_req(char *fname, char *commsyntax, hcge_cof *record, uvadd ncof, ha_cgeset *ha_set,uvdim nset);

int HaNDBBDMatOderPre(Mat A, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int nreg, int ntime, uvadd nvarele, uvadd *ha_eqadd,int *ha_rows,int *ha_cols, uvadd ndblock,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,ha_cgetype cntl6,PetscInt* ndbbdrank,PetscBool presol);
int HaDBBDMatOder(Mat A, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend, uvadd nvarele, uvadd *ha_eqadd,int *ha_rows,int *ha_cols, uvadd ndblock,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,ha_cgetype cntl6);
int HaDBBDMatOder1(Mat A, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend, uvadd nvarele, uvadd *ha_eqadd,int *ha_rows,int *ha_cols, uvadd ndblock,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,ha_cgetype cntl6);
int HaDBBDSol(Mat A, Vec b, ha_cgetype *x1, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, uvadd ndblock,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA, uvdim laD);
int HaDBBDParSol(Mat A, Vec b, ha_cgetype *x1, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, uvadd ndblock,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laD,PetscReal cntl3);//,bool iter
int HaNDBBDMatOder(Mat A, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int nreg, int ntime, uvadd nvarele, uvadd *ha_eqadd,int *ha_rows,int *ha_cols, uvadd ndblock,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,ha_cgetype cntl6,PetscInt* ndbbdrank,PetscBool presol);
int HaNDBBDParPre(Mat A, Vec b, ha_cgetype *x1, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, uvadd ndblock,uvadd nreg,uvadd ntime,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,uvdim laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol);//,bool iter
int HaNDBBDParSol(Mat A, Vec b, ha_cgetype *x1, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, uvadd ndblock,uvadd nreg,uvadd ntime,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,uvdim laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol);//,bool iter
int HaReduce(ha_cgetype *vecbivi,forint vecbivisize,PetscInt mpisize,PetscInt rank,PetscInt targetrank);
int HaReduceNoComp(ha_cgetype *vecbivi,forint vecbivisize,PetscInt mpisize,PetscInt rank,PetscInt targetrank);
int HaReducebivi(ha_cgetype *vecbivi,int *nz,int *biviindx1,int *bivinzrow,int *bivinzcol,int *bivinzrow0,int *bivinzcol0,int *bivinzrow1,int *bivinzcol1,int bivirowsize,int bivicolsize,PetscInt Dcolsize,PetscInt mpisize,PetscInt rank,PetscInt targetrank);
int HaSpaVecAdd(ha_cgetype *vecbivi,int *biviindx,int *nz0,ha_cgetype *vecbivi0,int *biviindx0,int j2);
int hcge_wdata(char *filename, char *newdatlogname, char *newdatfile,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele,hcge_cof *ha_cof,uvadd ncof,uvadd ncofele,hcge_cof *ha_var,uvadd nvar,uvadd nvarele, ha_cgevar *ha_cofvar);
//int HaMatGetSubMats(Mat A,uvadd indxbsize,PetscInt *indicesB,uvadd *ha_ndblocks,uvadd *begblock, uvadd *countvarintra1, uvadd *ha_cols,PetscInt mpisize,PetscInt rank,forint nmatin,PetscInt origrank,int **nzbs,int **bns,int **bms,int **irnbs,int **jcnbs,PetscScalar **valbs);

////Temp
//ha_cgetype ha_newfpcalshow(ha_cgevar *record,ha_cgeset *ha_set,ha_cgesetele *ha_setele,ha_cgesumele *ha_sumele,hcge_calvars *ha_calvar,int ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim, ha_cgetype zerodivide,int show);
uvadd hnew_biupd(PetscInt rank,char *fname,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele,ha_cgeexovar *ha_cgeshock,uvadd nvarele,int laA,uvdim subints,bool IsIni,int IsSplint,int nsteps);
bool spline(ha_cgetype* y,ha_cgetype* x,ha_cgetype sx0,ha_cgetype sxn,int size,ha_cgetype* w,int laA);

bool Johansen(PetscBool nohsl,PetscInt VecSize,Mat A,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,Mat B,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec vecb,Vec vece,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar **ha_cofvar2,uvadd ncofvar,uvadd ncofele,uvadd nvarele,ha_cgeexovar **ha_cgeshock2,uvadd alltimeset,uvadd allregset,uvadd nintraeq,uvdim matsol,PetscInt Istart,PetscInt Iend,  uvadd nreg, uvadd ntime, uvadd *ha_eqadd, uvadd ndblock, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,uvdim laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol,uvdim nesteddbbd,int localsize,PetscInt *ndbbddrank1,forint* indata,uvdim mc66,forint *ptx,struct timeval begintime,ha_cgetype **xcf2);
bool ModMidPoint(PetscBool nohsl,PetscInt VecSize,Mat* A,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,Mat* B,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec* vecb,Vec *vece,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar **ha_cofvar2,uvadd ncofvar,uvadd ncofele,uvadd nvarele,ha_cgeexovar **ha_cgeshock2,uvadd alltimeset,uvadd allregset,uvadd nintraeq,uvdim matsol,PetscInt Istart,PetscInt Iend,  uvadd nreg, uvadd ntime, uvadd *ha_eqadd, uvadd ndblock, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,uvdim laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol,uvdim nesteddbbd,int localsize,PetscInt *ndbbddrank1,forint* indata,uvdim mc66,forint *ptx,struct timeval begintime,uvdim subints,MPI_Fint fcomm,ha_cgetype **xcf2,int Isbiupd);
#endif // HA_CGEGLOBAL_H_INCLUDED

