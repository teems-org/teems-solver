#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

// I've replaced these with the ones from ha_cgeglobal.h to make sure it's ok
//#define CELLSIZE 20
//#define NAMESIZE 256//52 for dolder files
//#define MAXVARDIM 10 //maximum variable dimension
//#define HEADERSIZE 5
//#define TABREADLINE 20000//10000 for older files
//#define MAXSUPSET 12 // 10 for older files
//Single float forreal;

#define NAMESIZE 256
#define DATREADLINE 150000
#define TABREADLINE 20000//2536
#define HEADERSIZE 5
#define TABLINESIZE 20000//81
#define NOPERTINSUM 5
#define MAXVARDIM 10 //maximum variable dimension
#define MAXSUPSET 12
#define SORD 1
#define MAXSSIZE 187500000//1500000000/8

typedef int uvdim;
typedef long int uvadd;
typedef double forreal;

uvadd nvarele,nvar,nsetspace,nset;//,nset=142

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

// This was the old struct typedef that didn't work:
// typedef struct
// {
//   char cofname[NAMESIZE];
//   uvadd begadd;
//   uvdim size;
//   //uvdim dims[MAXVARDIM];//dimension size
//   //uvadd setbegadd;
//   uvadd setid[MAXVARDIM];
//   //uvadd setbegadd1[MAXVARDIM];
//   uvadd antidims[MAXVARDIM];
//   //char dimsets[MAXVARDIM][NAMESIZE];
//   uvadd matsize;
//   bool level_par;
//   bool change_real;
//   bool suplval;
// } hcge_cof ;

typedef float ha_floattype;
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


forreal* xc;//= (forreal *) calloc (1,sizeof(forreal));
ha_cgeset* ha_set;//= (ha_cgeset*) calloc (1,sizeof(ha_cgeset));
hcge_cof* ha_var;//= (hcge_cof) calloc (1,sizeof(hcge_cof));
ha_cgesetele* ha_setele;//= (ha_cgesetele*) calloc (1,sizeof(ha_cgesetele));
uvadd* varindxarray;//= (uvadd *) calloc (1,sizeof(uvadd));

int indx=0;
int strcmpchanged(const void* str1, const void* str2)
{
  return strcmp(*((char **)str1), *((char **)str2));
}

int main(int argc, char *argv[]) {
  long int i;
  uvdim j,l,dims[MAXVARDIM];
  size_t freadparm;
  double s;

  assert(argc == 2);

  char file[256]; 
  strcpy(file, argv[1]);

  FILE* data;
  for(i=0; i<256; i++)if(file[i]=='.' && file[i-1]=='l')break;
  file[i+1]='m';
  file[i+2]='d';
  file[i+3]='s';
  file[i+4]=0;
  uvadd modeldes[4];
  data = fopen(file, "rb");
  if ((data = fopen(file, "rb")) == NULL)
  {
    printf("Error opening .mds file\n");
    printf("Tried to open file with name: %s\n", file);
    exit(-1);
  }

  freadparm=fread(modeldes, sizeof(uvadd),4, data);
  fclose(data);
  nsetspace=modeldes[0];
  nvar=modeldes[1];
  nvarele=modeldes[2];
  nset=modeldes[3];

  printf("Contents of .mds file:\n");
  for (int j = 0; j < 4; j++) {
    printf("%ld\n", modeldes[j]);
  }

  xc=realloc(xc, sizeof(forreal)*nvarele);
  ha_set= realloc (ha_set,nset*sizeof(ha_cgeset));
  ha_setele= realloc (ha_setele,nsetspace*sizeof(ha_cgesetele));
  ha_var=realloc (ha_var,nvar*sizeof(hcge_cof));

  file[i+1]='b';
  file[i+2]='i';
  file[i+3]='n';
  file[i+4]=0;
  data = fopen(file, "rb");
  if ((data = fopen(file, "rb")) == NULL)
  {
    printf("Error opening .bin file\n");
  }
  freadparm=fread(xc, sizeof(forreal),nvarele, data);
  fclose(data);

  printf("Contents of .bin file:\n");
  for (int j = 0; j < nvarele; j++) {
    printf("%f\n", xc[j]);
  }

  file[i+1]='v';
  file[i+2]='a';
  file[i+3]='r';
  file[i+4]=0;
  // printf("file %s\n",file);
  data = fopen(file, "rb");
  if ((data = fopen(file, "rb")) == NULL)
  {
    printf("Error opening .var file\n");
  }
  freadparm=fread(ha_var, sizeof(hcge_cof),nvar, data);
  fclose(data);

  printf("Contents of .var file:\n");

  printf("%ld cofname:\n", nvar);
  for (int j = 0; j < nvar; j++) {
    printf("%s\n", ha_var[j].cofname);
  }

  printf("%ld begadd:\n", nvar);
  for (int j = 0; j < nvar; j++) {
    printf("%ld\n", ha_var[j].begadd);
  }

  printf("%ld size:\n", nvar);
  for (int j = 0; j < nvar; j++) {
    printf("%d\n", ha_var[j].size);
  }

  printf("%ld setid:\n", nvar);
  for (int j = 0; j < nvar; j++) {
    for (int k = 0; k < MAXVARDIM; k++) {
      printf("%ld ", ha_var[j].setid[k]); 
    }
    printf("\n");
  }

  printf("%ld antidims:\n", nvar);
  for (int j = 0; j < nvar; j++) {
    for (int k = 0; k < MAXVARDIM; k++) {
      printf("%ld ", ha_var[j].antidims[k]); 
    }
    printf("\n");
  }

  printf("%ld matsize:\n", nvar);
  for (int j = 0; j < nvar; j++) {
    printf("%ld\n", ha_var[j].matsize);
  }

  printf("%ld level_par:\n", nvar);
  for (int j = 0; j < nvar; j++) {
    printf("%d\n", ha_var[j].level_par);
  }

  printf("%ld change_real:\n", nvar);
  for (int j = 0; j < nvar; j++) {
    printf("%d\n", ha_var[j].change_real);
  }

  printf("%ld suplval:\n", nvar);
  for (int j = 0; j < nvar; j++) {
    printf("%d\n", ha_var[j].suplval);
  }

  printf("%ld gltype:\n", nvar);
  for (int j = 0; j < nvar; j++) {
    printf("%i\n", ha_var[j].gltype);
  }

  printf("%ld glval:\n", nvar);
  for (int j = 0; j < nvar; j++) {
    printf("%f\n", ha_var[j].glval);
  }

  file[i+1]='s';
  file[i+2]='e';
  file[i+3]='l';
  file[i+4]=0;
  data = fopen(file, "rb");
  if ((data = fopen(file, "rb")) == NULL)
  {
    printf("Error opening .sel file\n");
  }
  freadparm=fread(ha_setele, sizeof(ha_cgesetele),nsetspace, data);
  fclose(data);

  printf("Contents of .sel file:\n");

  printf("%ld setele:\n", nsetspace);
  for (int j = 0; j < nsetspace; j++) {
    printf("%s\n", ha_setele[j].setele);
  }

  printf("%ld setsh:\n", nsetspace);
  for (int j = 0; j < nsetspace; j++) {
    for (int k = 0; k < MAXVARDIM; k++) {
      printf("%d ", ha_setele[j].setsh[k]); 
    }
    printf("\n");
  }

  file[i+1]='s';
  file[i+2]='e';
  file[i+3]='t';
  file[i+4]=0;
  data = fopen(file, "rb");
  if ((data = fopen(file, "rb")) == NULL)
  {
    printf("Error opening .set file\n");
  }
  freadparm=fread(ha_set, sizeof(ha_cgeset),nset, data);
  fclose(data);
  printf("Contents of .set file:\n");

  printf("%ld header:\n", nset);
  for (int j = 0; j < nset; j++) {
    printf("%s\n", ha_set[j].header);
  }

  printf("%ld fileid:\n", nset);
  for (int j = 0; j < nset; j++) {
    printf("%i\n", ha_set[j].fileid);
  }

  printf("%ld setname:\n", nset);
  for (int j = 0; j < nset; j++) {
    printf("%s\n", ha_set[j].setname);
  }

  printf("%ld readele:\n", nset);
  for (int j = 0; j < nset; j++) {
    printf("%s\n", ha_set[j].readele);
  }

  printf("%ld begadd:\n", nset);
  for (int j = 0; j < nset; j++) {
    printf("%ld\n", ha_set[j].begadd);
  }

  printf("%ld size:\n", nset);
  for (int j = 0; j < nset; j++) {
    printf("%i\n", ha_set[j].size);
  }

  printf("%ld subsetid:\n", nset);
  for (int j = 0; j < nset; j++) {
    for (int k = 0; k < MAXSUPSET; k++) {
      printf("%d ", ha_set[j].subsetid[k]); 
    }
    printf("\n");
  }

  printf("%ld intertemp:\n", nset);
  for (int j = 0; j < nset; j++) {
    printf("%i\n", ha_set[j].intertemp);
  }

  printf("%ld intsup:\n", nset);
  for (int j = 0; j < nset; j++) {
    printf("%i\n", ha_set[j].intsup);
  }

  printf("%ld regional:\n", nset);
  for (int j = 0; j < nset; j++) {
    printf("%i\n", ha_set[j].regional);
  }

  printf("%ld regsup:\n", nset);
  for (int j = 0; j < nset; j++) {
    printf("%i\n", ha_set[j].regsup);
  }

}
