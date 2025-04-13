// Pham Van Ha & Tom Kompas
// Ensure that sum over time involves only intertemporal variables
// If otherwise switch back to normal MC66
// Currently no sub-reg set, can try different varreg to get opt: comment break in nestedmatvarread
// Currently logic fdim <> only in formula
#include <ha_cgeglobal.h>
//#include <ha_cgeiof.h>
//#include <ha_cgefparse.h>
//#include <ha_cgetab.h>
//#include <ha_cgename.h>
//#include "petscksp.h"
//struct point {
//forint nz,m,mpisize,nblock,nsbbdblocks;
//int *irn,*jcn;
//float *b,*values;
//};

static char help[] = "Solves a CGE model in parallel with KSP.\n\
           Input parameters include:\n\
           -None at the moment\n\n";

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args) {
  Vec      vecb,vece,x;  /* approx solution, RHS, exact solution */
  Mat      A,B;//,C,D;    /* linear system matrix */
  KSP      ksp;   /* linear solver context */
  PC       pc;   /* preconditionercontext */
  PetscRandom  rctx;   /* random number generator context */
  //PetscReal    norm;   /* norm of solution error */
  PetscInt  rank=0,mpisize,rank_hsl=0;//,wr_id,wr_color,wr_size
  PetscInt     VecSize=0,Istart=0,Iend=0,dnz=0,onz=0,dnzB=0,onzB=0,count,*onnz,*dnnz,*onnzB,*dnnzB,its,ndbbdrank=0;
  PetscErrorCode ierr;
  PetscBool   flg,presol;//,PreLoad = PETSC_FALSE;
  PetscScalar  value,zero=0;//,neg_one=-1;//*array;
  PetscLogDouble time0,time1;
  clock_t timestr,timeend,timemulti;
  struct timeval begintime, endtime,gettime_now;
  //struct timespec gettime_now,gettime_beg,gettime_end;
  //long int start_time=0;
  //double rep_time;
  size_t freadresult;
  uvadd lasize=1;
  uvadd i,j;
  uvadd j2=0,j1=0,j0=0,j3,j4,j5,j6;
  //VecScatter scatter; /* scatter context */
  //PetscViewer view_out;
  PetscInitialize(&argc,&args,(char *)0,help);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&mpisize);

  char processor_name[MPI_MAX_PROCESSOR_NAME+1];
  int name_len,name_len_max,name_beg,class_size,color,group_size,ha_id;
  MPI_Get_processor_name(processor_name, &name_len);
  printf("rank %d name len %d proc name %s\n",rank,name_len,processor_name);
  MPI_Allreduce(&name_len,&name_len_max,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
  name_len_max++;
  char *vec_pr_name=(char *) calloc (mpisize*name_len_max,sizeof(char));
  char *vec_pr_sname=(char *) calloc (mpisize*name_len_max,sizeof(char));
  name_beg=rank*name_len_max;
  for(i=name_beg; i<name_len+name_beg; i++)vec_pr_name[i]=processor_name[i-name_beg];
  vec_pr_name[i]='\0';
  //printf("rank %d name len %d proc name %s\n",rank,name_len,vec_pr_name+name_beg);
  MPI_Allreduce(vec_pr_name,vec_pr_sname,mpisize*name_len_max,MPI_CHAR,MPI_SUM,PETSC_COMM_WORLD);
  //if(rank==0)for(i=0;i<mpisize;i++)printf("snamessss %s\n",vec_pr_sname+i*name_len_max);
  for(i=0; i<name_len_max; i++)vec_pr_name[i]=vec_pr_sname[i];
  j=1;
  for(i=0; i<mpisize; i++) {
    for(j1=0; j1<j; j1++) {
      if(strncmp(vec_pr_sname+i*name_len_max,vec_pr_name+j1*name_len_max,name_len_max)==0) {
        break;
      }
    }
    if(j1==j) {
      //printf("rank %d sname %s name %s\n",rank,vec_pr_sname+i*name_len_max,vec_pr_name+j1*name_len_max);
      for(j6=0; j6<name_len_max; j6++)vec_pr_name[j*name_len_max+j6]=vec_pr_sname[j1*name_len_max+j6];
      j++;
    }
  }
  class_size=j;
  for(i=0; i<class_size; i++) {
    if(strncmp(processor_name,vec_pr_name+i*name_len_max,name_len_max)==0) {
      color=i;
      break;
    }
  }
  free(vec_pr_name);
  free(vec_pr_sname);
  MPI_Comm_split(PETSC_COMM_WORLD,color,rank,&HA_COMM);
  MPI_Comm_rank( HA_COMM, &ha_id);
  MPI_Comm_size(HA_COMM,&group_size);
  if(ha_id==group_size-1)color=1;
  else color=0;
  MPI_Comm_split(PETSC_COMM_WORLD,color,rank,&HA1_COMM);

  //MPI_Comm WR_COMM;
  //MPI_Comm_split_type(PETSC_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,&WR_COMM);
  //MPI_Comm_rank( WR_COMM, &wr_id);
  //ierr = PetscGetCPUTime(&time0);
  //CHKERRQ(ierr);
  gettimeofday(&begintime, NULL);
  bool sbbd_overrid=false;
  PetscBool sbbd_overuser=false,nohsl=false;//,iter=false;
  char ch='y';
  if(rank==0) {
    printf("***********\nWarnings:\n***********\nA. This software was written by Pham Van Ha and Tom Kompas for their own research. Use it at your own risk!\nB. The shocks to only 2 dimensions variable are different from GEMPACK. In fact the matrix needs to be transposed.\nIntertemporal variables (equation) must be declared with minimum dimension to minimise\nnet cut (e.g. qo(COM,REG,TIME) must be capital(REG,TIME)=qo(\"capital\",REG,TIME)\nFor NDBBD, LA in Di is not impotant! The procedure automatically determine rank!\nTake care with text file created in windows, dos environment when using gedit. Gedit will use cr for new line! Bug!\nLaA, LaDi are very important. They control the size of temporary files. Enter smallest posible!!!\nDo you agree with the term (y/n)?\n");
    //ch = getchar();
  }
  //if(ch!='y') return 0;
  MPI_Barrier(PETSC_COMM_WORLD);
  //**************************************************************************************
  //****************************** READ SET ELEMENT***************************************
  //**************************************************************************************
  char tabfile[TABREADLINE],newtabfile[TABREADLINE]="_temp_tab_file",newtabfile1[TABREADLINE]="_temp_tab_new_file",closure[TABREADLINE],shock[TABREADLINE],filename[TABREADLINE],longname[TABREADLINE],vname[NAMESIZE],copyline[TABREADLINE],regset[NAMESIZE];
  char tempfilenam[255],tempchar[255],solmed[NAMESIZE],solchar[255];
  int niodata=0,nj,mem_fac=0,noutdata=0,nsoldata=0,nowrites=0;
  uvadd nsetspace=0,dcount,ndblock=0,netcut=0,ndblock1,nreg=0,ntime=0;
  uvdim nset=0,vsize,dim1,nlength=0,matsol=0,laA=2,laDi=2,laD=2,nsbbdblocks=2,nesteddbbd=0,mc66=0,subints=1,subindx,StoIter=1;
  uvadd alltimeset=-1,allregset=-1;
  PetscReal cntl6=0,cntl3;
  if(rank<10) {
    strcat(newtabfile,"000");
    strcat(newtabfile1,"000");
  }
  if(rank<100&&rank>9) {
    strcat(newtabfile,"00");
    strcat(newtabfile1,"00");
  }
  if(rank<1000&&rank>99) {
    strcat(newtabfile,"0");
    strcat(newtabfile1,"0");
  }
  sprintf(tempchar, "%d",rank);
  //printf("temchar %s\n",tempchar);
  //printf("sof petscscalar %d sof double %d\n",sizeof(PetscScalar),sizeof(double));
  strcat(newtabfile,tempchar);
  strcat(newtabfile1,tempchar);
  strcat(newtabfile,".tab");
  strcat(newtabfile1,".tab");
  //printf("Enter the command file name: ");
  //scanf("%s",&filename);
  //printf("fname %s\n",filename);
  PetscOptionsGetInt(NULL,NULL,"-matsol",&matsol,NULL);//0 MA48 (nproc must be 1) 1 SBBD 2 DBBD, if >=1 should enable reg or time. First reg set will be regorder in var. Note INCL(4) in hsl_mp48ss
  if(matsol==2)nohsl=true;
  if(matsol==3)nohsl=true;
  isLinux=0;
  PetscOptionsGetInt(NULL,NULL,"-laA",&laA,NULL);
  if(laA==0)laA=2;
  PetscOptionsGetInt(NULL,NULL,"-laD",&laD,NULL);
  if(laD==0)laD=2;
  PetscOptionsGetInt(NULL,NULL,"-laDi",&laDi,NULL);
  if(laDi==0)laDi=2;
  PetscOptionsGetInt(NULL,NULL,"-withmc66",&mc66,NULL);
  //printf("mc66 %d\n",mc66);
  PetscOptionsGetInt(NULL,NULL,"-step1",&step1,NULL);
  if(step1==0)step1=2;
  PetscOptionsGetInt(NULL,NULL,"-step2",&step2,NULL);
  if(step2==0)step2=4;
  PetscOptionsGetInt(NULL,NULL,"-step3",&step3,NULL);
  if(step3==0)step3=8;
  kindx1=step1/(double)2;
  i=(uvadd)step1/2;
  if(kindx1!=i){
    //odd
    kindx1=step2/(double)2;
    i=(uvadd)step2/2;
    if(kindx1==i){
      printf("Error!!! steps must be all odd or even\n");
      return 0;
    }
    kindx1=step3/(double)2;
    i=(uvadd)step3/2;
    if(kindx1==i){
      printf("Error!!! steps must be all odd or even\n");
      return 0;
    }
  }else{
    //even
    kindx1=step2/(double)2;
    i=(uvadd)step2/2;
    if(kindx1!=i){
      printf("Error!!! steps must be all odd or even\n");
      return 0;
    }
    kindx1=step3/(double)2;
    i=(uvadd)step3/2;
    if(kindx1!=i){
      printf("Error!!! steps must be all odd or even\n");
      return 0;
    }
  }
  
  kindx1=step2/(double)step1;
  step2=(PetscInt)step2/step1;
  kindx2=step3/(double)step1;
  step3=(PetscInt)step3/step1;
//   if(kindx1!=step2){
//     printf("Error!!! step2 must be a multiple of step1\n");
//     return 0;
//   }
//   if(kindx1!=step2){
//     printf("Error!!! step3 must be a multiple of step1\n");
//     return 0;
//   }
  smallthreads=0;
  medthreads=0;
  mymaxnumthrd=1;
  PetscOptionsGetInt(NULL,NULL,"-maxthreads",&mymaxnumthrd,NULL);
  if(mymaxnumthrd>1&&mymaxnumthrd<=omp_get_max_threads( )){
    omp_set_num_threads(mymaxnumthrd);
    //printf("Please press Enter!!!\n");
    //getchar();
  }else{
    printf("Max Threads Num = %d\nPlease set OMP_NUM_THREADS <=maxthreads!!! I am setting it to 1!\n",mymaxnumthrd);
    mymaxnumthrd=1;
    omp_set_num_threads(mymaxnumthrd);
  }
  PetscOptionsGetInt(NULL,NULL,"-medthreads",&smallthreads,NULL);
  PetscOptionsGetInt(NULL,NULL,"-smllthreads",&smallthreads,NULL);
  if(smallthreads==0)smallthreads=mymaxnumthrd;
  if(medthreads==0)medthreads=mymaxnumthrd;
  PetscOptionsGetInt(NULL,NULL,"-nsubints",&subints,NULL);
  PetscOptionsGetInt(NULL,NULL,"-nsbbdblocks",&nsbbdblocks,NULL);
  PetscOptionsGetInt(NULL,NULL,"-nesteddbbd",&nesteddbbd,NULL);
  PetscOptionsGetInt(NULL,NULL,"-nowrites",&nowrites,NULL);
  PetscOptionsGetBool(NULL,NULL,"-presol",&presol,NULL); /* preparation for next solution */
  PetscOptionsGetReal(NULL,NULL,"-cntl_6",&cntl6,NULL); /* CNTL6 in Mat Order */
  PetscOptionsGetReal(NULL,NULL,"-cntl_3",&cntl3,NULL);/*Iterative threshold */
  PetscOptionsGetInt(NULL,NULL,"-ndbbd_bl_rank",&ndbbdrank,NULL);/*Override default rank for last block in NDBBD method. Read text file. >0 text column. Use with care*/
  PetscOptionsGetInt(NULL,NULL,"-isLinux",&isLinux,NULL);/*For read and fread, slightly better performance for linux*/
  PetscOptionsGetInt(NULL,NULL,"-stoiter",&StoIter,NULL);
  isLinux=0;//permanently set to zero, errors in large system!!!!
  PetscInt *ndbbddrank1= (PetscInt *) calloc (ndbbdrank,sizeof(PetscInt));
  PetscOptionsGetString(NULL,NULL,"-nestfile",filename,TABREADLINE,&flg);
  if (!flg) {
    strcpy(filename,"./ndbbd_drank.csv");//orani03.cmf");
  }
  //printf("matsol %s\n",filename);
  if (flg)intreadCSV(filename,ndbbddrank1,ndbbdrank);//
  //for(i=0;i<ndbbdrank;i++)printf("i %d rank %d\n",i,ndbbddrank1[i]);
  printf("matsol %d\n",matsol);
  PetscOptionsGetString(NULL,NULL,"-cmdfile",filename,TABREADLINE,&flg);
  if (!flg) {
    strcpy(filename,"./reg.cmf");//orani03.cmf");
  }
  PetscOptionsGetBool(PETSC_NULL,NULL, "-enable_time", &sbbd_overuser,PETSC_NULL);/* Overrid MC66 ordering */
  //PetscOptionsGetBool(PETSC_NULL,NULL, "-enable_iter", &iter,PETSC_NULL);/* Overrid MC66 ordering */
  //PetscOptionsGetBool(PETSC_NULL, "-no_hsl", &nohsl,PETSC_NULL);/* No HSL Petsc solution */
  regset[0]='\0';
  PetscOptionsGetString(NULL,NULL,"-regset",regset,NAMESIZE,&flg);
  //printf("matsol1 %s\n",regset);
  if(regset[0]!='\0')for(i=0; i<NAMESIZE; i++) {
      regset[i]=tolower((int)regset[i]);
    }
  PetscOptionsGetString(NULL,NULL,"-solmed",solmed,NAMESIZE,&flg);
  //printf("matsol %d\n",matsol);
  if (!flg) {
    strcpy(solmed,"Mmid");//orani03.cmf");
  }
  int solmethod;
  if(strcmp(solmed,"Mmid")==0)solmethod=1;
  if(strcmp(solmed,"Johansen")==0)solmethod=10;
  if(strcmp(solmed,"Stochastic")==0)solmethod=20;
  if(strcmp(solmed,"StoSim")==0)solmethod=21;
  if(strcmp(solmed,"NoSol")==0)solmethod=100;
  printf("Sol med %d regset %s\n",solmethod,regset);

  #pragma omp parallel private(i)
  {
  i=0;
  }

  if(nohsl) {
    rank_hsl=rank;
  }
  else {
    rank_hsl=0;
  }
  if(rank==rank_hsl) {
    mem_fac=1;
  }
  /*FILE * filehandle;
  char linetry[TABREADLINE];
  filehandle = fopen("mo.cmf","r");
  fgets(linetry,TABREADLINE,filehandle);
  fclose(filehandle);
  ierr = PetscGetCPUTime(&time1);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Calculation of variables time %f\n",time1-time);
  ierr = PetscFinalize();
  CHKERRQ(ierr);
  return 0;*/
  char *readitem=NULL;
  if(rank==0) {
    printf("OK!!!\n");
    niodata=hcge_niodata(filename,"iodata");
    if(niodata==-1)return 0;
    printf("OK!!!\n");
    noutdata=hcge_niodata(filename,"outdata");
    nsoldata=hcge_niodata(filename,"soldata");
  }
  if(nohsl) {
    MPI_Bcast(&niodata,sizeof(int), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&noutdata,sizeof(int), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nsoldata,sizeof(int), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  hcge_iodata *iodata= (hcge_iodata *) calloc (niodata+noutdata+nsoldata,sizeof(hcge_iodata));
  //printf("rank %d OKK1\n",rank);
  if(rank==rank_hsl) {
    //printf("rank %d nio %d OK!!\n",rank,niodata);
    //if(rank==0) {
    hcge_rcmd(filename,niodata,iodata,tabfile,closure,shock);
    //printf("rank1 %d nio %d OK!!\n",rank,niodata);
    for (nj=0; nj<niodata+noutdata+nsoldata; nj++) printf("rank %d logname %s fname %s\n",rank,iodata[nj].logname,iodata[nj].filname);
    //printf("f %s tab %s cls %s shf %s\n",filename,tabfile,closure,shock);
    //printf("OKgh!\n");
    if(hcge_wtab(tabfile,newtabfile)==-1)return 0;
    printf("OK1!\n");
  }
  printf("rank %d OKK1\n",rank);

  strcpy(tabfile,newtabfile);
  if(rank==0)nset=ha_cgenset(tabfile);
  //}
  if(nohsl) {
    MPI_Bcast(iodata,niodata*sizeof(hcge_iodata), MPI_BYTE,0, PETSC_COMM_WORLD);
    //MPI_Bcast(tabfile,TABREADLINE*sizeof(char), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(closure,TABREADLINE*sizeof(char), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(shock,TABREADLINE*sizeof(char), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nset,sizeof(uvdim), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  ha_cgeset *ha_set= (ha_cgeset *) calloc (nset,sizeof(ha_cgeset));
  for(i=0; i<nset; i++) {
    ha_set[i].subsetid[0]=i;
    for(j=0; j<MAXSUPSET; j++)ha_set[i].subsetid[j]=-1;
  }
  //printf("nset %d\n",nset);
  if(rank==0) {
    ha_cgerset(tabfile,niodata,iodata, ha_set,nset);
    //printf("OKK!!!\n");
    hcge_rinterset(tabfile,niodata,iodata, ha_set,nset);
    //for (j=0; j<nset; j++) printf("setname %s header %s setsize %d readele %s inter %d\n",ha_set[j].setname,ha_set[j].header,ha_set[j].size,ha_set[j].readele,ha_set[j].intertemp);
    for (i=0; i<nset; i++) {
      ha_set[i].begadd=nsetspace;
      nsetspace=nsetspace+ha_set[i].size;
    }
    if(regset[0]!='\0')for (i=0; i<nset; i++)if(strcmp(regset,ha_set[i].setname)==0) {
          allregset=i;
          break;
        }
  }
  if(nohsl) {
    MPI_Bcast(ha_set,nset*sizeof(ha_cgeset), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&allregset,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nsetspace,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  printf("rank %d regset %s indx %ld\n",rank,regset,allregset);
  ha_cgesetele *ha_setele= (ha_cgesetele *) calloc (nsetspace,sizeof(ha_cgesetele));
  printf("nset %d nsetspace %ld\n",nset,nsetspace);
  //dcount=0;
  for (i=0; i<nsetspace; i++)for (j=0; j<MAXSUPSET; j++)ha_setele[i].setsh[j]=-1;
  if(rank==0) {
    for (i=0; i<nset; i++) {
      strcpy(vname,ha_set[i].header);
      nlength=0;
      while (vname[nlength] != '\0') {
        nlength++;
      }
      //printf("i %ld nl %d vname %s\n",i,nlength,ha_set[i].header);
      if (nlength>0) {
        ha_cgerdvar1(vname,iodata[ha_set[i].fileid].filname,&vsize,longname,&dim1);
        //printf("dim1 %d\n",dim1);
        ha_cgemvar1 *matvar1= (ha_cgemvar1 *) calloc (dim1,sizeof(ha_cgemvar1));
        ha_cgermvar1(vname,iodata[ha_set[i].fileid].filname,dim1,matvar1);
        //printf("OK!!!\n");
        for (j=0; j<dim1; j++) {
          //printf("matva %s len %d dim1 %d j %d\n",matvar1[j].ch,strlen(matvar1[j].ch),dim1,j);
          nj=0;
          while(matvar1[j].ch[nj]!='\0') {
            matvar1[j].ch[nj]=tolower((int) matvar1[j].ch[nj]);
            nj++;
          }
          strncpy(ha_setele[j+ha_set[i].begadd].setele,matvar1[j].ch,strlen(matvar1[j].ch));
          ha_setele[j+ha_set[i].begadd].setsh[0]=j;
          //printf("name %s setsh %d\n",ha_set[i].setname,j);
          //ha_set[i].subsetid[0]=i;
          //strcpy(ha_set[i].supersetname[0],ha_set[i].setname);
          //printf("setele %s setsh %d\n",ha_setele[j+dcount].setele,ha_setele[j+dcount].setsh);
        }
        //dcount=dcount+dim1;
        free(matvar1);
      }
      else {
        if (ha_set[i].readele[0]=='-'&&ha_set[i].readele[1]==',') {
          ha_setminus(ha_setele, ha_set,nset,i);
          //dcount=dcount+ha_set[i].size;
        }
        else {
          if (ha_set[i].readele[0]=='+'&&ha_set[i].readele[1]==',') {
            ha_setplus(ha_setele, ha_set,nset,i);
            //dcount=dcount+ha_set[i].size;
          }
          else {
            if (ha_set[i].readele[0]=='^'&&ha_set[i].readele[1]==',') {
              ha_setunion(ha_setele, ha_set,nset,i);
              //dcount=dcount+ha_set[i].size;
            }
            else {
              if(ha_set[i].readele[0]=='=') {
                dim1=ha_set[i].size;
                strcpy(copyline,ha_set[i].readele);
                //strcat(copyline,"=");
                //printf("i %d line %s dim1 %d\n",i,copyline,dim1);
                while (ha_cgefrstr(copyline," ", ""));
                readitem = strtok(copyline,"=");
                //printf("j11111 %d %s\n",j1,readitem);
                readitem = strtok(readitem,"=");
                //printf("j11111 %d %s\n",j1,readitem);
                j1=atoi(readitem);
                //printf("j11111 %d %s\n",j1,readitem);
                for (j=0; j<dim1; j++) {
                  strcpy(ha_setele[j+ha_set[i].begadd].setele,ha_setele[j+ha_set[j1].begadd].setele);
                  ha_setele[j+ha_set[i].begadd].setsh[0]=j;
                  //printf("setele %s\n",ha_setele[j+ha_set[j1].begadd].setele);
                }
              }
              else {
                dim1=ha_set[i].size;
                strcpy(copyline,ha_set[i].readele);
                strcat(copyline,",");
                //printf("i %d line %s\n",i,copyline);
                while (ha_cgefrstr(copyline," ", ""));
                readitem = strtok(copyline,",");
                //printf("line1 %s dim1 %d dcount %d\n",readitem,dim1,ha_set[i].begadd);
                strcpy(ha_setele[ha_set[i].begadd].setele,readitem);
                ha_setele[ha_set[i].begadd].setsh[0]=0;
                for (j=1; j<dim1; j++) {
                  readitem = strtok(NULL,",");
                  strcpy(ha_setele[j+ha_set[i].begadd].setele,readitem);
                  ha_setele[j+ha_set[i].begadd].setsh[0]=j;
                  //printf("name %s setsh %d\n",ha_set[i].setname,j);
                  //ha_set[i].subsetid[0]=i;
                  //strcpy(ha_set[i].supersetname[0],ha_set[i].setname);
                }
                //dcount=dcount+dim1;
              }
            }
          }
        }
      }

    }
    //for (i=0;i<nset;i++) for (j=0;j<ha_set[i].size;j++) printf("i %d setname %s begad %d size %d setele %s short %d sup1 %d sup2 %d sub3 %d\n",i,ha_set[i].setname,ha_set[i].begadd,ha_set[i].size,ha_setele[ha_set[i].begadd+j].setele,ha_setele[ha_set[i].begadd+j].setsh[0],ha_setele[ha_set[i].begadd+j].setsh[1],ha_setele[ha_set[i].begadd+j].setsh[2],ha_set[i].subsetid[2]);
    //printf("OK!!!\n");
    ha_cgersubset(tabfile, ha_setele, ha_set,nset);
    j2=1;
    while(j2==1)for(i=1; i<MAXSUPSET; i++)ha_cgesubsetchck(ha_setele,ha_set,nset,&j2); //printf("check %d\n",i);}
    if(sbbd_overuser) {
      alltimeset=ha_cgeralltime(ha_set,nset);
    }
    //for (i=0;i<nset;i++)for (j2=0;j2<MAXSUPSET;j2++)printf("i %ld name %s size %d subset %d\n",i,ha_set[i].setname,ha_set[i].size,ha_set[i].subsetid[j2]);
    //for (i=0;i<nset;i++) for (j=0;j<ha_set[i].size;j++) printf("i %d setname %s begad %d size %d setele %s short %d sup1 %d sup2 %d sub3 %d\n",i,ha_set[i].setname,ha_set[i].begadd,ha_set[i].size,ha_setele[ha_set[i].begadd+j].setele,ha_setele[ha_set[i].begadd+j].setsh[0],ha_setele[ha_set[i].begadd+j].setsh[1],ha_setele[ha_set[i].begadd+j].setsh[3],ha_set[i].subsetid[3]);
    //alltimeset=-1;
    if(alltimeset>=0||allregset>=0) {
      if(alltimeset>=0&&allregset>=0) {
        ndblock=ha_set[alltimeset].size*ha_set[allregset].size;
        nreg=ha_set[allregset].size;
        ntime=ha_set[alltimeset].size;
      }
      if(alltimeset>=0&&allregset>=0&&nesteddbbd==1) {
        ndblock=ha_set[alltimeset].size*(ha_set[allregset].size+1);
        nreg=ha_set[allregset].size;
        ntime=ha_set[alltimeset].size;
      }
      if(alltimeset<0&&allregset>=0) {
        ndblock=ha_set[allregset].size;
        nreg=ndblock;
      }
      if(alltimeset>=0&&allregset<0) {
        ndblock=ha_set[alltimeset].size;
        ntime=ndblock;
      }
    }
    if(allregset>=0) {
      ha_set[allregset].regional=true;
      //for(j=1;j<MAXSUPSET;j++)if(ha_set[allregset].subsetid[j]>-1)ha_set[ha_set[allregset].subsetid[j]].regional=true;else break;
      for(i=0; i<nset; i++) {
        for(j=1; j<MAXSUPSET; j++) {
          if(ha_set[i].subsetid[j]>-1) {
            if(ha_set[ha_set[i].subsetid[j]].regional) {
              ha_set[i].regional=true;
              ha_set[i].regsup=j;
              break;
            }
          }
          else break;
        }
      }
      //}
    }
    if(alltimeset>=0) {
      for(i=0; i<nset; i++) {
        for(j=1; j<MAXSUPSET; j++) {
          if(ha_set[i].subsetid[j]>-1) {
            if(ha_set[ha_set[i].subsetid[j]].intertemp) {
              ha_set[i].intsup=j;
              break;
            }
          }
          else break;
        }
      }
    }
    ndblock1=ndblock;
    printf("rank %d alltime set %ld allreg %ld size %ld\n",rank,alltimeset,allregset,ndblock);
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  if(nohsl) {
    MPI_Bcast(&alltimeset,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&allregset,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(ha_set,nset*sizeof(ha_cgeset), MPI_BYTE,0, PETSC_COMM_WORLD);
    //for(i=0;i<nsetspace;i++)MPI_Bcast(ha_setele+i,sizeof(ha_cgesetele), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(ha_setele,nsetspace*sizeof(ha_cgesetele), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&ndblock1,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&ntime,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nreg,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  if(ntime!=ndbbdrank&&nesteddbbd==1) {
    printf("Please set ndbbdrank = ntime %ld!\n",ntime);
    return 0;
  }
  ndblock=ndblock1;
  printf("rank %d ndblock %ld allreg %ld\n",rank,ndblock,allregset);
  //for (i=0;i<nset;i++) for (j=0;j<ha_set[i].size;j++) printf("setname %s id %d begad %d size %d setele %s short %d supid1 %d sup1 %d supid2 %d sub2 %d supid3 %d sub3 %d\n",ha_set[i].setname,i,ha_set[i].begadd,ha_set[i].size,ha_setele[ha_set[i].begadd+j].setele,ha_setele[ha_set[i].begadd+j].setsh[0],ha_set[i].subsetid[1],ha_setele[ha_set[i].begadd+j].setsh[1],ha_set[i].subsetid[2],ha_setele[ha_set[i].begadd+j].setsh[2],ha_set[i].subsetid[3],ha_setele[ha_set[i].begadd+j].setsh[3]);
  //return 0;
  //for (j=0;j<nset;j++) printf("setname %s sub %d\n",ha_set[j].setname,ha_set[j].subsetid);

  //**************************************************************************************
  //****************************** END READ SET ELEMENT***********************************
  //**************************************************************************************

  //**************************************************************************************
  //****************************** READ COEFFICIENT NAME**********************************
  //**************************************************************************************
  char commsyntax[NAMESIZE];
  strcpy(commsyntax,"coefficient");
  uvadd ncof=0,ncofele=0,ncof1,ncofele1;
  if(rank==0) {
    printf("tabfile %s\n",tabfile);
    ncof=ha_cgencof(tabfile,commsyntax);
    printf("tabfile %s ncof %ld\n",tabfile,ncof);
    ncof1=ncof;
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  printf("rank %d ncof %ld\n",rank,ncof);
  if(nohsl)MPI_Bcast(&ncof1,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
  ncof=ncof1;
  //return 0;
  //printf("ncof %d\n",ncof);
//  ha_cgecof *ha_cof= (ha_cgecof *) calloc (ncof,sizeof(ha_cgecof));//recycle ha_cgeset
  hcge_cof *ha_cof= (hcge_cof *) calloc (ncof,sizeof(hcge_cof));//recycle ha_cgeset
  //printf("nset %d\n",nset);
  if(rank==0) {
    ncofele=hcge_rcof(tabfile,commsyntax,ha_cof,ncof,ha_set,nset);
    if(ncofele==-1)return 0;
    ncofele1=ncofele;
  }
  if(nohsl) {
    MPI_Bcast(&ncofele1,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
    //for(i=0;i<ncof;i++)MPI_Bcast(ha_cof+i,sizeof(hcge_cof), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(ha_cof,ncof*sizeof(hcge_cof), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  ncofele=ncofele1;
  printf("rank %d ncofele %ld\n",rank,ncofele);
  //for (i=0;i<ncof;i++) printf("size %d begadd %d cname %s d1 %d d2 %d d3 %d d4 %d d5 %d d1n %s d2n %s d3n %s d4n %s d5n %s\n",ha_cof[i].size,ha_cof[i].begadd,ha_cof[i].cofname,ha_cof[i].dims[0],ha_cof[i].dims[1],ha_cof[i].dims[2],ha_cof[i].dims[3],ha_cof[i].dims[4],ha_cof[i].dimsets[0],ha_cof[i].dimsets[1],ha_cof[i].dimsets[2],ha_cof[i].dimsets[3],ha_cof[i].dimsets[4]);

  //for (i=0;i<ncof;i++) printf("size %d begadd %d cname %s d1 %d d2 %d d3 %d d4 %d d1n %s d2n %s d3n %s d4n %s\n",ha_cof[i].size,ha_cof[i].begadd,ha_cof[i].cofname,ha_cof[i].d1,ha_cof[i].d2,ha_cof[i].d3,ha_cof[i].d4,ha_cof[i].d1name,ha_cof[i].d2name,ha_cof[i].d3name,ha_cof[i].d4name);
  //for (i=0;i<ncofele;i++) printf("cofele %s\n",ha_cofele[i].cofname);
  //return 0;
  //**************************************************************************************
  //****************************** END READ COEFFICIENT NAME******************************
  //**************************************************************************************

  //**************************************************************************************
  //****************************** READ VARIABLE NAME*************************************
  //**************************************************************************************
  strcpy(commsyntax,"variable");
  uvadd nvar=0,nvarele=0,nvar1,nvarele1;
  if(rank==0) {
    nvar=ha_cgencof(tabfile,commsyntax);
    nvar1=nvar;
  }
  if(nohsl)MPI_Bcast(&nvar1,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
  nvar=nvar1;
//  ha_cgecof *ha_var= (ha_cgecof *) calloc (nvar,sizeof(ha_cgecof));//recycle ha_cgeset
  hcge_cof *ha_var= (hcge_cof *) calloc (nvar,sizeof(hcge_cof));//recycle ha_cgeset
  bool *var_inter= (bool *) calloc (nvar,sizeof(bool));//recycle ha_cgeset
  printf("nvarele %ld\n",nvarele);
  if(rank==0) {
    hcge_defvar(tabfile,ha_var,nvar);
    nvarele=hcge_rvar(tabfile,commsyntax,ha_var,nvar,ha_set,nset);
    nvarele1=nvarele;
  }
  printf("nvarele %ld\n",nvarele);
  if(nohsl) {
    MPI_Bcast(&nvarele1,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
    //for(i=0;i<nvar;i++)MPI_Bcast(ha_var+i,sizeof(hcge_cof), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(ha_var,nvar*sizeof(hcge_cof), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  nvarele=nvarele1;
  if(rank==0){
  for (i=0;i<nvar;i++){
    //printf("i %ld %s\n",i,ha_var[i].cofname);
    if(ha_var[i].cofname[1]=='_'){
    for (j=0;j<ncof;j++)if (strcmp(ha_cof[j].cofname,ha_var[i].cofname+2)==0){
      printf("Error!!! Same variable and coefficient names are not supported in this version, even with p_ or c_\nPlease change the name of coefficient %s or variable %s\n",ha_cof[j].cofname,ha_var[i].cofname);
      return 0;
    }
    }
  }
  }
//for (i=0;i<nvar;i++){printf("var %s ",ha_var[i].cofname);for (j=0;j<ha_var[i].size;j++) {printf("dim %d %d ",j,ha_var[i].setid[j]);}printf("\n");}
  printf("nvarele %ld\n",nvarele);
  ha_cgevar *ha_cofvar= (ha_cgevar *) calloc ((ncofele+nvarele),sizeof(ha_cgevar));
  ha_cgecofele *ha_cofele= (ha_cgecofele *) calloc (ncofele,sizeof(ha_cgecofele));
  ha_cgecofele *ha_varele= (ha_cgecofele *) calloc (nvarele,sizeof(ha_cgecofele));
  printf("rankasd %d nvar %ld\n",rank,nvar);
  //if(rank==0) {
  //j=0;
  //for (i=0; i<nvar; i++) {
  //  j+=ha_var[i].matsize;
  //}
  //printf("j %d varele %d\n",j,nvarele);
  //for (i=0;i<nvar;i++) printf("size %d\n",ha_var[i].size);
  //for (i=0;i<nvar;i++) printf("size %d begadd %d cname %s d1 %d d2 %d d3 %d d4 %d d5 %d d1n %s d2n %s d3n %s d4n %s d5n %s anti1 %d anti2 %d matsize %d\n",ha_var[i].size,ha_var[i].begadd,ha_var[i].cofname,ha_var[i].dims[0],ha_var[i].dims[1],ha_var[i].dims[2],ha_var[i].dims[3],ha_var[i].dims[4],ha_var[i].dimsets[0],ha_var[i].dimsets[1],ha_var[i].dimsets[2],ha_var[i].dimsets[3],ha_var[i].dimsets[4],ha_var[i].antidims[0],ha_var[i].antidims[1],ha_var[i].matsize);
  //printf("nvarele %d\n",nvarele);
  if(rank==0) {
    hcge_rcofele(ha_cof,ncof,ha_set,nset,ha_cofele);
  }
  //if(nohsl)for(i=0;i<ncofele;i++)MPI_Bcast(ha_cofele+i,sizeof(ha_cgecofele), MPI_BYTE,0, PETSC_COMM_WORLD);
  if(nohsl) {
    if(ncofele*sizeof(ha_cgecofele)>1500000000) {
      j1=1500000000/sizeof(ha_cgecofele);
      i=ncofele/j1;
      for(j=0; j<i; j++) {
        MPI_Bcast(ha_cofele+j*j1,j1*sizeof(ha_cgecofele), MPI_BYTE,0, PETSC_COMM_WORLD);
      }
      i=ncofele-j*j1;
      MPI_Bcast(ha_cofele+j*j1,i*sizeof(ha_cgecofele), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(ha_cofele,ncofele*sizeof(ha_cgecofele), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
    //MPI_Bcast(ha_cofele,ncofele*sizeof(ha_cgecofele), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  printf("rank %d ncofele %ld\n",rank,ncofele);

  if(rank==0)hcge_rcofele(ha_var,nvar,ha_set,nset,ha_varele);
  printf("rank %d OK!!!\n",rank);
  //for (i=0;i<nvarele;i++) printf("varele %s\n",ha_varele[i].cofname);

  if(rank==rank_hsl)hcge_wvar(tabfile,newtabfile1,ha_var,nvar);
  //printf("OK!!!\n");
  strcpy(tabfile,newtabfile1);
  //}
  if(nohsl) {
    //MPI_Bcast(tabfile,TABREADLINE*sizeof(char), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(ha_set,nset*sizeof(ha_cgeset), MPI_BYTE,0, PETSC_COMM_WORLD);
    //MPI_Bcast(ha_var,nvar*sizeof(hcge_cof), MPI_BYTE,0, PETSC_COMM_WORLD);
    //MPI_Bcast(ha_varele,nvarele*sizeof(ha_cgecofele), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  //printf("tabfile %s\n",tabfile);
  //strcpy(tabfile,"_temp_tab_new_file1.tab");
  //**************************************************************************************
  //****************************** END READ VARIABLE NAME*********************************
  //**************************************************************************************

  //**************************************************************************************
  //********************* READ VARIABLE, COEFFICIENT VALUE FROM FILE**********************
  //**************************************************************************************
  if(rank==0) {
    strcpy(commsyntax,"read");
    //long int nread;
    //nread=
    if(hcge_readff(tabfile,niodata,iodata,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_cofele,ncofele,ha_var,nvar,ha_varele,nvarele)==-1)return 0;
  }
  printf("rank %d OK???\n",rank);
//  FILE *ofpt;
//  ofpt = fopen("xvart.txt", "w");
//  j0=0;
//  for(i=0; i<ncof; i++) {
//    for(j=0;j<ha_cof[i].matsize;j++){
//        fprintf(ofpt, "Cof Name: %s Cof Val %f\n",ha_cof[i].cofname,ha_cofele[j0].cofval);
//        j0++;
//    }
//  }
//  fclose(ofpt);
  //for (i=0;i<nvar;i++)for (j=0;j<ha_var[i].matsize;j++) printf("varele %s val %lf\n",ha_var[i].cofname,ha_varele[ha_var[i].begadd+j].cofval);
  //for (i=0;i<ncof;i++)for (j=0;j<ha_cof[i].matsize;j++) printf("varele %s val %lf\n",ha_cof[i].cofname,ha_cofele[ha_var[i].begadd+j].cofval);
  //for (i=0;i<ncofele;i++) printf("cofele %s val %lf\n",ha_cofele[i].cofname,ha_cofele[i].cofval);
  //printf("nread %d\n",nread);
  //return 0;
  //**************************************************************************************
  //********************* END READ VARIABLE, COEFFICIENT VALUE FROM FILE******************
  //**************************************************************************************

  //**************************************************************************************
  //***************** CALCULATE VARIABLE, COEFFICIENT VALUE FROM FORMULA******************
  //**************************************************************************************

  uvadd nexo=0,nexo1;//,nshock,nform,
  //printf("ncofvar %ld\n",ncofele+nvarele);
  ha_cgevar *ha_cofvar1;
  //printf("ncofvar %ld size of ele %d\n",ncofele+nvarele,sizeof(ha_cgevar));
  //ha_cgevar *ha_cofvar1= (ha_cgevar *) calloc ((ncofele+nvarele),sizeof(ha_cgevar));
  if(rank==0) {
    for (i=0; i<ncofele; i++) {
      //strcpy(ha_cofvar[i].varname,ha_cofele[i].cofname);
      ha_cofvar[i].varval=ha_cofele[i].cofval;
    }
    //printf("ncofvar %ld\n",ncofele+nvarele);
    for (i=ncofele; i<nvarele+ncofele; i++) {
      //strcpy(ha_cofvar[i].varname,ha_varele[i-ncofele].cofname);
      ha_cofvar[i].varval=ha_varele[i-ncofele].cofval;
    }
  }
  printf("rank %d ncofvar %ld\n",rank,ncofele+nvarele);
  free(ha_cofele);
  //for (j=0; j<nvarele; j++) {
  //  strcpy(ha_cgeshock[j].ExoName,ha_varele[j].cofname);
  //}
  free(ha_varele);
  ha_cgeexovar *ha_cgeshock= (ha_cgeexovar *) calloc (nvarele,sizeof(ha_cgeexovar));
  if(rank==0) {
    strcpy(commsyntax,"exogenous");
    //printf("OK!!!!\n");
    nexo=hcge_rexo(closure,commsyntax,ha_cgeshock,ha_var,nvar,ha_set,nset,ha_setele);
    nexo1=nexo;
    //printf("nexo %d\n",nexo);
    strcpy(commsyntax,"shock");
    //return 0;
    //nshock=
    if(hcge_rshock(shock,commsyntax,ha_cgeshock,nvarele,ha_var,nvar,ha_set,nset,ha_setele,subints)==-1)return 0;
    //printf("nexo %d nvarele %d\n",nexo,nvarele);
    //for (i=0; i<nvar; i++)for (j=0; j<ha_var[i].matsize; j++) printf("varele %s exo %d exoindx %ld val %lf\n",ha_var[i].cofname,ha_cgeshock[ha_var[i].begadd+j].ShockId,ha_cgeshock[ha_var[i].begadd+j].ExoIndx,ha_cgeshock[ha_var[i].begadd+j].ShockVal);
    //for (i=0; i<nvarele; i++) printf("%ld val %ld\n",i,ha_cgeshock[i].ExoIndx);
  }
  printf("rank %d OK???\n",rank);
  if(nohsl) {
    //for(i=0;i<nvarele;i++)MPI_Bcast(ha_cgeshock+i,sizeof(ha_cgeexovar), MPI_BYTE,0, PETSC_COMM_WORLD);
    //printf("rank %d size %ld\n",rank,nvarele*sizeof(ha_cgeexovar));
    if(nvarele*sizeof(ha_cgeexovar)>1500000000) {
      j1=1500000000/sizeof(ha_cgeexovar);
      i=nvarele/j1;
      for(j=0; j<i; j++) {
        //printf("rank %d j %d size %ld send %ld nvarele %ld\n",rank,j,j1,(j+1)*j1,nvarele);
        MPI_Bcast(ha_cgeshock+j*j1,j1*sizeof(ha_cgeexovar), MPI_BYTE,0, PETSC_COMM_WORLD);
      }
      i=nvarele-j*j1;
      //printf("rank %d j %d size %ld send %ld\n",rank,j,j1,i);
      MPI_Bcast(ha_cgeshock+j*j1,i*sizeof(ha_cgeexovar), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(ha_cgeshock,nvarele*sizeof(ha_cgeexovar), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
    //MPI_Bcast(ha_cgeshock,nvarele*sizeof(ha_cgeexovar), MPI_BYTE,0, PETSC_COMM_WORLD);
    MPI_Bcast(&nexo1,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  nexo=nexo1;
  strcpy(commsyntax,"formula");
  bool IsIni=true;
  printf("OK???\n");
  //return 1;
  //for (i=23; i<24; i++) printf("%d cofvar %s val %lf\n",i-ncofele,ha_cofvar[i].varname,ha_cofvar[i].varval);
  //nform=
  if(rank==0) {
    hnew_calcff(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,IsIni);
  }
  //ierr = PetscGetCPUTime(&time1);
  //CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"Calculation of variables time %f\n",time1-time0);
  //CHKERRQ(ierr);
  gettimeofday(&endtime, NULL);
  if(rank==0)printf("Calculation of variables time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
  //ierr = PetscGetCPUTime(&time0);
  //CHKERRQ(ierr);
  if(nohsl) { //Overcome MPI_Bcast limit
    if((nvarele+ncofele)*sizeof(ha_cgevar)>1500000000) {
      j1=1500000000/sizeof(ha_cgevar);
      i=(nvarele+ncofele)/j1;
      for(j=0; j<i; j++) {
        MPI_Bcast(ha_cofvar+j*j1,j1*sizeof(ha_cgevar), MPI_BYTE,0, PETSC_COMM_WORLD);
      }
      i=nvarele+ncofele-j*j1;
      MPI_Bcast(ha_cofvar+j*j1,i*sizeof(ha_cgevar), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(ha_cofvar,(nvarele+ncofele)*sizeof(ha_cgevar), MPI_BYTE,0, PETSC_COMM_WORLD);
    }
    //MPI_Barrier(PETSC_COMM_WORLD);
  }
  //printf("rank!!!!!!!!!!!!!!!!! %d cover %f nohsl %d send size %ld var size %d\n",rank,ha_cofvar[16895628].varval,nohsl,(nvarele+ncofele),sizeof(ha_cgevar));
  //return 0;
  //for(i=0;i<ncof;i++)if(strcmp(ha_cof[i].cofname,"fc1")==0)for(j4=0;j4<ha_cof[i].matsize;j4++)printf("vname %s %.8f\n",ha_cof[i].cofname,ha_cofvar[ha_cof[i].begadd+j4].varval);
  //for(i=0;i<ncof;i++)if(strcmp(ha_cof[i].cofname,"vendwwld")==0)for(j4=0;j4<ha_cof[i].matsize;j4++)printf("vname %s %f\n",ha_cof[i].cofname,ha_cofvar[ha_cof[i].begadd+j4].varval);
  //for (i=0;i<ncof;i++)for (j=0;j<ha_cof[i].matsize;j++) printf("cofele %s val %.8f\n",ha_cof[i].cofname,ha_cofvar[ha_cof[i].begadd+j].varval);
  //for (i=0;i<nvar;i++)for (j=0;j<ha_var[i].matsize;j++) printf("varele %s val %lf\n",ha_var[i].cofname,ha_cofvar[ncofele+ha_var[i].begadd+j].varval);
  //for (i=0; i<nvarele; i++) ha_cofvar1[i].varval=ha_cofvar[i].varval;
//  FILE *ofpt;
//  ofpt = fopen("xvart.txt", "w");
//  j0=0;
//  for(i=0; i<ncof; i++) {
//    for(j=0;j<ha_cof[i].matsize;j++){
//        fprintf(ofpt, "Cof Name: %s Cof Val %f\n",ha_cof[i].cofname,ha_cofvar[j0].varval);
//        j0++;
//    }
//  }
//  fclose(ofpt);
//  FILE *ofpt;
//  ofpt = fopen("xvart.txt", "w");
//  for(i=0;i<ncofele;i++){
//    fprintf(ofpt, "Cof Name: %s Cof Val %f\n",ha_cofvar[i].varname,ha_cofvar[i].varval);
//  }
//  for(i=0;i<nvarele;i++){
//    if(ha_cgeshock[i].ShockId==0){fprintf(ofpt, "Var Name: %s Var Val %f\n",ha_cofvar[ncofele+i].varname,ha_cofvar[i+ncofele].varval);}
//    else {fprintf(ofpt, "Var Name: %s Var Val %f\n",ha_cofvar[ncofele+i].varname,ha_cgeshock[i].ShockVal);}
//  }
//  fclose(ofpt);
  //printf("nform %d\n",nform);
  //for (i=0; i<(nvarele+ncofele); i++) printf("%ld cofvar %lf\n",i-ncofele,ha_cofvar[i].varval);
  //return 0;
  //**************************************************************************************
  //****************** END CALCULATE VARIABLE, COEFFICIENT VALUE FROM FORMULA*************
  //**************************************************************************************
  //ierr = PetscGetCPUTime(&time1);
  //CHKERRQ(ierr);
  gettimeofday(&begintime, NULL);
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"Broadcast of variables time %f\n",time1-time0);
  if(rank==0)printf("Broadcast of variables time %f\n",(begintime.tv_sec - endtime.tv_sec)+((double)(begintime.tv_usec - endtime.tv_usec))/ 1000000);
  //CHKERRQ(ierr);
  //ierr = PetscGetCPUTime(&time0);
  //CHKERRQ(ierr);
  //return 0;
  //**************************************************************************************
  //****************************** MATRIX FROM FORMULA************************************
  //**************************************************************************************
  VecSize = (PetscInt) nvarele-nexo;
  //if(rank==rank_hsl) {
  PetscPrintf(PETSC_COMM_SELF,"VecSize %d exo %d\n",VecSize,nexo);
  //}
  //PetscPrintf(PETSC_COMM_SELF,"Ist!\n");
  //printf("OK1!!!\n");
  //ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,VecSize,&x);
  //CHKERRQ(ierr);
  //printf("OK1!!!\n");
  //ierr = VecSetFromOptions(x);
  //CHKERRQ(ierr);
  //ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,VecSize,&e);//nexo
  //CHKERRQ(ierr);
  //ierr = VecSetFromOptions(e);
  //CHKERRQ(ierr);
  strcpy(commsyntax,"equation");
  //PetscPrintf(PETSC_COMM_SELF,"here\n");
  //printf("OK1!!!\n");
  //PetscSynchronizedPrintf(PETSC_COMM_SELF,"Istart: %d Iend %d\n",Istart,Iend);
  uvadd neq=0,neq1;
  if(rank==0) {
    neq=ha_cgencof(tabfile,commsyntax);
    neq1=neq;
  }
  if(nohsl) {
    MPI_Bcast(&neq1,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
  }
  neq=neq1;
  if(rank==rank_hsl) {
    printf("neq %ld\n",neq);
  }
  uvadd *countvarintra1= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
  hcge_cof *ha_eq= (hcge_cof *) calloc (neq,sizeof(hcge_cof));//recycle ha_cgeset
  bool *ha_eqint= (bool *) calloc (neq,sizeof(bool));//recycle ha_cgeset
  uvdim *orderintra= (uvdim *) malloc (nvar*sizeof(uvdim));
  uvdim *orderreg= (uvdim *) malloc (nvar*sizeof(uvdim));
  for(i=0; i<nvar; i++) {
    orderintra[i]=-1;
    orderreg[i]=-1;
  }
  uvdim *ha_eqtime= (uvdim *) malloc (neq*sizeof(uvdim));//recycle ha_cgeset
  uvdim *ha_eqreg= (uvdim *) malloc (neq*sizeof(uvdim));//recycle ha_cgeset
  for(i=0; i<neq; i++) {
    ha_eqtime[i]=-1;
    ha_eqreg[i]=-1;
  }
  uvadd nintraendovar,summat;
  printf("rank %d\n",rank);
  if(rank==rank_hsl) {
    //for(i=0;i<neq;i++)printf("i %d neq %d\n",i,ha_eqint[i]);
    if(nesteddbbd==1)NestedMatvarRead(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,var_inter,ha_eq,ha_eqint,ha_eqtime,ha_eqreg,allregset,alltimeset,orderintra,orderreg);
    else NewMatvarRead(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,var_inter,ha_eq,ha_eqint,ha_eqtime,ha_eqreg,allregset,alltimeset,orderintra,orderreg);
    if(alltimeset>=0||allregset>=0)for(i=0; i<neq; i++)ha_eqint[i]=!ha_eqint[i];
    //for(i=0;i<neq;i++)printf("rank %d i %ld eqname %s neq %d eqreg %d\n",rank,i,ha_eq[i].cofname,ha_eqint[i],ha_eqreg[i]);
    //for(i=0;i<nvar;i++)printf("rank %d var %s int %d size %d\n",rank,ha_var[i].cofname,var_inter[i],ha_var[i].matsize);
    //for(i=0;i<nvar;i++)printf("rank %d i %d var %s ordintra %d ordreg %d varinter %d\n",rank,i,ha_var[i].cofname,orderintra[i],orderreg[i],var_inter[i]);
    //summat=0;
    //for(i=0;i<nvar;i++)if(var_inter[i])summat+=ha_var[i].matsize;
    //printf("var inter size %d\n",summat);
  }
  //for(i=0;i<nvar;i++)printf("%s var inter size %d odrint %d ordreg %d\n",ha_var[i].cofname,var_inter[i],orderintra[i],orderreg[i]);
  switch (nesteddbbd) {
  case 1 :
    if(!(alltimeset>=0&&allregset>=0))printf("Not a intertemporal regional CGE model!\n");
    uvadd *countvarintra= (uvadd *) calloc (ndblock,sizeof(uvadd));
    j3=0;
    //j5=0;
    for (i=0; i<nvar; i++) {
      for (j=0; j<ha_var[i].matsize; j++) {
        if(!ha_cgeshock[j3+j].ShockId) {
          if(!var_inter[i]) {
            j0=j;
            j2=-1;
            for(j1=0; j1<orderintra[i]+1; j1++) {
              j2=j0/ha_var[i].antidims[j1];
              j0-=j2*ha_var[i].antidims[j1];
            }
            j0=j;
            j4=-1;
            for(j1=0; j1<orderreg[i]+1; j1++) {
              j4=j0/ha_var[i].antidims[j1];
              j0-=j4*ha_var[i].antidims[j1];
            }
            //if(orderintra[i]>-1){
            if(j4>-1)if(ha_set[ha_var[i].setid[orderreg[i]]].regsup>0)j4=ha_setele[ha_set[ha_var[i].setid[orderreg[i]]].begadd+j4].setsh[ha_set[ha_var[i].setid[orderreg[i]]].regsup];
            if(ha_set[ha_var[i].setid[orderintra[i]]].intsup>0)j2=ha_setele[ha_set[ha_var[i].setid[orderintra[i]]].begadd+j2].setsh[ha_set[ha_var[i].setid[orderintra[i]]].intsup];
            if(orderreg[i]>-1)countvarintra[j2*(nreg+1)+j4]++;
            else countvarintra[j2*(nreg+1)+nreg]++;
            //}
          }
        }//else j5++;
      }
      //printf("rank %d i %ld hamat %ld j3 %ld j5 %ld\n",rank,i,ha_var[i].matsize,j3,j5);
      j3+=ha_var[i].matsize;
    }
    if(rank==rank_hsl) {
      j2=countvarintra[0];
      countvarintra[0]=0;
      nintraendovar=j2;
    }
    for (i=1; i<ndblock; i++) {
      nintraendovar+=countvarintra[i];
      j3=countvarintra[i];
      countvarintra[i]=countvarintra[i-1]+j2;
      countvarintra1[i]=countvarintra[i];
      //printf("rank %d count1 %ld\n",rank,countvarintra[i]);
      j2=j3;
    }
    j0=0;
    j1=0;
    j2=0;
    j3=0;
    j4=nintraendovar;
    for (i=0; i<nvar; i++) {
      for (j=0; j<ha_var[i].matsize; j++) {
        j5=j3+j;
        if(!ha_cgeshock[j5].ShockId) {
          if(!var_inter[i]) {
            j0=j;
            j2=-1;
            for(j1=0; j1<orderintra[i]+1; j1++) {
              j2=j0/ha_var[i].antidims[j1];
              j0-=j2*ha_var[i].antidims[j1];
            }
            j0=j;
            j6=-1;
            for(j1=0; j1<orderreg[i]+1; j1++) {
              j6=j0/ha_var[i].antidims[j1];
              j0-=j6*ha_var[i].antidims[j1];
            }
            if(ha_set[ha_var[i].setid[orderintra[i]]].intsup>0)j2=ha_setele[ha_set[ha_var[i].setid[orderintra[i]]].begadd+j2].setsh[ha_set[ha_var[i].setid[orderintra[i]]].intsup];
            if(j6>-1)if(ha_set[ha_var[i].setid[orderreg[i]]].regsup>0)j6=ha_setele[ha_set[ha_var[i].setid[orderreg[i]]].begadd+j6].setsh[ha_set[ha_var[i].setid[orderreg[i]]].regsup];
            if(orderreg[i]>-1) {
              ha_cgeshock[j5].ExoIndx=countvarintra[j2*(nreg+1)+j6];
              countvarintra[j2*(nreg+1)+j6]++;
            }
            else {
              ha_cgeshock[j5].ExoIndx=countvarintra[j2*(nreg+1)+nreg];
              countvarintra[j2*(nreg+1)+nreg]++;
            }
          }
          else {
            ha_cgeshock[j5].ExoIndx=j4;
            //if(ha_cgeshock[j5].ExoIndx==246026)printf("OKj5!!!!!!!!!!!!!!! %s\n",ha_var[i].cofname);
            j4++;
          }
        }
      }
      j3+=ha_var[i].matsize;
    }
    //for(i=0;i<ndblock;i++)printf("counintra %d\n",countvarintra[i]);
    if(rank==rank_hsl) {
      countvarintra1[ndblock]=countvarintra[ndblock-1];
      //printf("count1 %ld\n",countvarintra1[ndblock]);
    }
    j1=0;
    for (i=0; i<nvarele; i++) {
      if (ha_cgeshock[i].ShockId) {
        ha_cgeshock[i].ExoIndx+=j1;
        j1++;
      }
    }
    //for(i=0;i<ndblock+1;i++)printf("counintra %d counintra1 %d\n",countvarintra[i],countvarintra1[i]);
    free(countvarintra);

    break;
  default :
    if(alltimeset>=0) {
      /*for (i=0; i<nvar; i++) {
        j3=0;
        j4=0;
        for(j=0; j<ha_var[i].size; j++) {
          for(j0=1; j0<nset; j0++) {
            if(strcmp(ha_set[ha_var[i].setid[j]].setname,ha_set[j0].setname)==0) {
              if(ha_set[j0].intertemp) {
                j3++;
                orderintra[i]=j;
              }
              break;
            }
          }
          if(allregset>=0)if(strcmp(ha_set[ha_var[i].setid[j]].setname,regset)==0) {
              orderreg[i]=j;
              j4++;
            }
        }
        if(j4==0&&allregset>=0) {
          var_inter[i]=true;
        }
        if(j3==0||ha_var[i].size==0) {
          var_inter[i]=true;
        }
      }*/
      //for (i=0; i<nvar; i++) if(var_inter[i])printf("intervar %s\n",ha_var[i].cofname);
      uvadd *countvarintra= (uvadd *) calloc (ndblock,sizeof(uvadd));
      j3=0;
      for (i=0; i<nvar; i++) {
        for (j=0; j<ha_var[i].matsize; j++) {
          if(!ha_cgeshock[j3+j].ShockId) {
            if(!var_inter[i]) {
              j0=j;
              for(j1=0; j1<orderintra[i]+1; j1++) {
                j2=j0/ha_var[i].antidims[j1];
                j0-=j2*ha_var[i].antidims[j1];
              }
              if(allregset>=0) {
                j0=j;
                for(j1=0; j1<orderreg[i]+1; j1++) {
                  j4=j0/ha_var[i].antidims[j1];
                  j0-=j4*ha_var[i].antidims[j1];
                }
                if(ha_set[ha_var[i].setid[orderreg[i]]].regsup>0)j4=ha_setele[ha_set[ha_var[i].setid[orderreg[i]]].begadd+j4].setsh[ha_set[ha_var[i].setid[orderreg[i]]].regsup];
                countvarintra[j2*nreg+j4]++;
              }
              else {
                countvarintra[j2]++;
              }
            }
          }
        }
        j3+=ha_var[i].matsize;
      }
      //for (i=0; i<ha_set[alltimeset].size; i++) printf("intrasize %d\n",countvarintra[i]);
      if(rank==rank_hsl) {
        j2=countvarintra[0];
        countvarintra[0]=0;
        nintraendovar=j2;
      }
      for (i=1; i<ndblock; i++) {
        nintraendovar+=countvarintra[i];
        j3=countvarintra[i];
        countvarintra[i]=countvarintra[i-1]+j2;
        countvarintra1[i]=countvarintra[i];
        j2=j3;
      }
      //for (i=0; i<ha_set[alltimeset].size; i++) if(!var_inter[i])printf("intrasize %d total %d\n",countvarintra[i],nintraendovar);
      j0=0;
      j1=0;
      j2=0;
      j3=0;
      j4=nintraendovar;
      //for (i=0; i<ha_set[alltimeset].size; i++) printf("intrasize %d\n",countvarintra[i]);
      //printf("nintraendovar %d\n",nintraendovar);
      for (i=0; i<nvar; i++) {
        //if(var_inter[i])printf("nam %s matsize %d\n",ha_var[i].cofname,ha_var[i].matsize);
        for (j=0; j<ha_var[i].matsize; j++) {
          j5=j3+j;
          if(!ha_cgeshock[j5].ShockId) {
            if(!var_inter[i]) {
              j0=j;
              for(j1=0; j1<orderintra[i]+1; j1++) {
                j2=j0/ha_var[i].antidims[j1];
                j0-=j2*ha_var[i].antidims[j1];
              }
              if(allregset>=0) {
                j0=j;
                for(j1=0; j1<orderreg[i]+1; j1++) {
                  j6=j0/ha_var[i].antidims[j1];
                  j0-=j6*ha_var[i].antidims[j1];
                }
                if(ha_set[ha_var[i].setid[orderreg[i]]].regsup>0)j6=ha_setele[ha_set[ha_var[i].setid[orderreg[i]]].begadd+j6].setsh[ha_set[ha_var[i].setid[orderreg[i]]].regsup];
                ha_cgeshock[j5].ExoIndx=countvarintra[j2*ha_set[allregset].size+j6];
                countvarintra[j2*ha_set[allregset].size+j6]++;
              }
              else {
                ha_cgeshock[j5].ExoIndx=countvarintra[j2];
                countvarintra[j2]++;
              }
            }
            else {
              //if(j4==nintraendovar)printf("nam %s j5-1 %d j5-2 %d\n",ha_var[i].cofname,ha_cgeshock[j5-1].ExoIndx,ha_cgeshock[j5-2].ExoIndx);
              ha_cgeshock[j5].ExoIndx=j4;
              j4++;
            }
          }
        }
        j3+=ha_var[i].matsize;
      }
      //printf("j4 %d\n",j4-nintraendovar);
      if(rank==rank_hsl) {
        countvarintra1[ndblock]=countvarintra[ndblock-1];
      }
      j1=0;
      for (i=0; i<nvarele; i++) {
        if (ha_cgeshock[i].ShockId) {
          ha_cgeshock[i].ExoIndx+=j1;
          j1++;
        }
        //if (i>361799&&i<362400) printf("indx %d\n",ha_cgeshock[i].ExoIndx);
      }
      //for(i=0;i<ndblock;i++)printf("counintra %d\n",countvarintra[i]);
      free(countvarintra);
    }
    else {
      if(allregset>=0) {
        /*uvdim *orderreg= (uvdim *) calloc (nvar,sizeof(uvdim));
        for (i=0; i<nvar; i++) {
          j4=0;
          for(j=0; j<ha_var[i].size; j++) if(strcmp(ha_set[ha_var[i].setid[j]].setname,regset)==0) {
              orderreg[i]=j;
              j4++;
              break;
            }
          if(j4==0) {
            var_inter[i]=true;
          }
        }*/
        //for (i=0; i<nvar; i++) if(var_inter[i])printf("intervar %s\n",ha_var[i].cofname);
        //for (i=0; i<nvar; i++) if(!var_inter[i])printf("intravar %s order %d\n",ha_var[i].cofname,orderreg[i]);
        uvadd *countvarintra= (uvadd *) calloc (ndblock,sizeof(uvadd));
        j3=0;
        for (i=0; i<nvar; i++) {
          for (j=0; j<ha_var[i].matsize; j++) {
            if(!ha_cgeshock[j3+j].ShockId) {
              if(!var_inter[i]) {
                //printf("var %s\n",ha_var[i].cofname);
                j0=j;
                for(j1=0; j1<orderreg[i]+1; j1++) {
                  j4=j0/ha_var[i].antidims[j1];
                  j0-=j4*ha_var[i].antidims[j1];
                }
                if(ha_set[ha_var[i].setid[orderreg[i]]].regsup>0)j4=ha_setele[ha_set[ha_var[i].setid[orderreg[i]]].begadd+j4].setsh[ha_set[ha_var[i].setid[orderreg[i]]].regsup];
                countvarintra[j4]++;
              }
            }
          }
          j3+=ha_var[i].matsize;
        }
        //for (i=0; i<ha_set[allregset].size; i++) printf("intrasize %d\n",countvarintra[i]);
        if(rank==rank_hsl) {
          j2=countvarintra[0];
          countvarintra[0]=0;
          nintraendovar=j2;
        }
        for (i=1; i<ndblock; i++) {
          nintraendovar+=countvarintra[i];
          j3=countvarintra[i];
          countvarintra[i]=countvarintra[i-1]+j2;
          countvarintra1[i]=countvarintra[i];
          j2=j3;
        }
        //for (i=0; i<ha_set[alltimeset].size; i++) if(!var_inter[i])printf("intrasize %d total %d\n",countvarintra[i],nintraendovar);
        j0=0;
        j1=0;
        j2=0;
        j3=0;
        j4=nintraendovar;
        //for (i=0; i<ha_set[alltimeset].size; i++) printf("intrasize %d\n",countvarintra[i]);
        //printf("nintraendovar %d\n",nintraendovar);
        for (i=0; i<nvar; i++) {
          //if(var_inter[i])printf("nam %s matsize %d\n",ha_var[i].cofname,ha_var[i].matsize);
          for (j=0; j<ha_var[i].matsize; j++) {
            j5=j3+j;
            if(!ha_cgeshock[j5].ShockId) {
              if(!var_inter[i]) {
                j0=j;
                for(j1=0; j1<orderreg[i]+1; j1++) {
                  j6=j0/ha_var[i].antidims[j1];
                  j0-=j6*ha_var[i].antidims[j1];
                }
                if(ha_set[ha_var[i].setid[orderreg[i]]].regsup>0)j6=ha_setele[ha_set[ha_var[i].setid[orderreg[i]]].begadd+j6].setsh[ha_set[ha_var[i].setid[orderreg[i]]].regsup];
                ha_cgeshock[j5].ExoIndx=countvarintra[j6];
                countvarintra[j6]++;
              }
              else {
                //if(j4==nintraendovar)printf("nam %s j5-1 %d j5-2 %d\n",ha_var[i].cofname,ha_cgeshock[j5-1].ExoIndx,ha_cgeshock[j5-2].ExoIndx);
                ha_cgeshock[j5].ExoIndx=j4;
                j4++;
              }
            }
          }
          j3+=ha_var[i].matsize;
        }
        //printf("j4 %d\n",j4-nintraendovar);
        if(rank==rank_hsl) {
          countvarintra1[ndblock]=countvarintra[ndblock-1];
        }
        j1=0;
        for (i=0; i<nvarele; i++) {
          if (ha_cgeshock[i].ShockId) {
            ha_cgeshock[i].ExoIndx+=j1;
            j1++;
          }
          //if (i>361799&&i<362400) printf("indx %d\n",ha_cgeshock[i].ExoIndx);
        }
        free(countvarintra);
      }
      else {
        j1=0;
        j2=0;
        for (i=0; i<nvar; i++) {
          for (j=0; j<ha_var[i].matsize; j++) {
            j3=j0+j;

            if (!ha_cgeshock[j3].ShockId) {
              ha_cgeshock[j3].ExoIndx+=j2;
              j2++;
            }
            if (ha_cgeshock[j3].ShockId) {
              ha_cgeshock[j3].ExoIndx+=j1;
              j1++;
            }
          }
          j0+=ha_var[i].matsize;
          //if (i>361799&&i<362400) printf("indx %d\n",ha_cgeshock[i].ExoIndx);
        }
      }
    }
  }
  free(orderintra);
  free(orderreg);
  //for (i=0; i<nvar; i++)for (j=0; j<ha_var[i].matsize; j++) printf("varele %s exo %d exoindx %ld val %lf\n",ha_var[i].cofname,ha_cgeshock[ha_var[i].begadd+j].ShockId,ha_cgeshock[ha_var[i].begadd+j].ExoIndx,ha_cgeshock[ha_var[i].begadd+j].ShockVal);
  //printf("shock %s OK!!!\n",shock);
  //for (i=0; i<nvar; i++) if(var_inter[i])printf("intervar %s\n",ha_var[i].cofname);

  free(var_inter);
  strcpy(commsyntax,"equation");
  uvadd *ha_eqadd= (uvadd *) calloc (VecSize,sizeof(uvadd));//recycle ha_cgeset
  //uvadd *ha_eqtimematsize= (uvadd *) calloc (neq,sizeof(uvadd));//recycle ha_cgeset
  uvadd *ha_eqtimesbegad= (uvadd *) calloc (neq,sizeof(uvadd));
  uvadd *ha_eqregsbegad= (uvadd *) calloc (neq,sizeof(uvadd));
  uvadd *counteq= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
  uvadd *counteqnoadd= (uvadd *) calloc (ndblock,sizeof(uvadd));
  uvadd nintraeq=0;
  //return 0;
  //printf("t %d r %d\n",alltimeset,allregset);
  if(alltimeset>=0&&allregset<0) {
    //printf("t %d r %d\n",alltimeset,allregset);
    //if(rank==rank_hsl) {
    //hcge_req(tabfile,commsyntax,ha_eq,neq,ha_set,nset);
    //}
    //for(i=0; i<neq; i++) {
    //  printf("eq name: %s\n",ha_eq[i].cofname);
    //  for(j=0; j<ha_eq[i].size+1; j++)printf("eq dim %d %s\n",j,ha_eq[i].dimsets[j]);
    //}
    //for(j0=0; j0<neq; j0++)for(j=0; j<ha_eq[j0].size; j++)printf("interset %s\n",ha_set[ha_eq[j0].setid[j]].setname);
    for(i=0; i<neq; i++) {
      j3=1;
      //j1=0;
      //printf("eq size %d\n",ha_eq[i].size);
      /*for(j=0; j<ha_eq[i].size; j++)for(j0=0; j0<nset; j0++)if(strcmp(ha_set[ha_eq[i].setid[j]].setname,ha_set[j0].setname)==0) {
            if(ha_set[j0].intertemp) {
              ha_eqtime[i]=j;
              ha_eqtimesbegad[i]=ha_set[j0].begadd;
              j1++;
            }
            //printf("j %d interset %d\n",j,ha_set[j0].intertemp);
            ha_eq[i].setid[j]=j0;
            j3*=ha_set[j0].size;
            break;
          }
      if(j1>0) {
        ha_eqint[i]=true;
      }
      ha_eq[i].matsize=j3;*/
      if(ha_eqtime[i]>-1)ha_eqtimesbegad[i]=ha_set[ha_eq[i].setid[ha_eqtime[i]]].begadd;
      //printf("timebegadd %d tbegadd %d\n",ha_set[ha_eq[i].setid[ha_eqtime[i]]].begadd,ha_set[alltimeset].begadd);
      /*for(j=0; j<ha_eq[i].size; j++){
            j3*=ha_set[ha_eq[i].setid[j]].size;
          }
      ha_eq[i].matsize=j3;*/
      //printf("matsize %d\n",ha_eq[i].matsize);
      //ha_eqtimematsize[i]=ha_eq[i].matsize;
      ha_eq[i].antidims[ha_eq[i].size-1]=1;
      if(ha_eq[i].size>1) {
        for (j2=ha_eq[i].size-2; j2>-1; j2--) {
          ha_eq[i].antidims[j2]=ha_eq[i].antidims[j2+1]*ha_set[ha_eq[i].setid[j2+1]].size;
        }
      }
    }
    //j1=0;
    //for(i=0; i<neq; i++) j1+=ha_eq[i].matsize;
    //printf("Number of equations %d\n",j1);
    j3=0;
    for (i=0; i<neq; i++) {
      if(ha_eqint[i]) {
        for (j=0; j<ha_eq[i].matsize; j++) {
          j0=j;
          for(j1=0; j1<ha_eqtime[i]+1; j1++) {
            j2=j0/ha_eq[i].antidims[j1];
            j0-=j2*ha_eq[i].antidims[j1];
          }
          if(ha_eqtime[i]>-1)if(ha_eq[i].setid[ha_eqtime[i]]==alltimeset)counteq[ha_setele[ha_eqtimesbegad[i]+j2].setsh[0]]++;
            else {
              for(j3=1; j3<MAXSUPSET; j3++)if(ha_set[ha_eq[i].setid[ha_eqtime[i]]].subsetid[j3]=alltimeset)break;
              counteq[ha_setele[ha_eqtimesbegad[i]+j2].setsh[j3]]++;
            }
        }
      }//else printf("i %d eq %s\n",i,ha_eq[i].cofname);
      //j3+=ha_eq[i].matsize;
    }
    //for (i=0; i<ha_set[alltimeset].size; i++)printf("intraeq %d\n",counteq[i]);
    if(rank==rank_hsl) {
      j2=counteq[0];
      counteqnoadd[0]=counteq[0];
      counteq[0]=0;
      nintraeq=j2;
    }
    for (i=1; i<ndblock; i++) {
      counteqnoadd[i]=counteq[i];
      nintraeq+=counteq[i];
      j3=counteq[i];
      counteq[i]=counteq[i-1]+j2;
      j2=j3;
    }
    if(rank==rank_hsl) {
      counteq[ndblock]=VecSize;
      netcut=VecSize-countvarintra1[ndblock];
      //for (i=0; i<ndblock+1; i++) if(ha_eqint[i])printf("intraeq1 %d total %d\n",counteq[i],nintraeq);
    }
  }
  //if(rank==rank_hsl) for (i=0; i<ndblock+1; i++) printf("c %d c1 %d netcut %d\n",counteq[i],countvarintra1[i],netcut);

  if(alltimeset<0&&allregset>=0) {
    //if(rank==rank_hsl) {
    //hcge_req(tabfile,commsyntax,ha_eq,neq,ha_set,nset);
    //}
    //for(i=0; i<neq; i++) {
    //  printf("eq name: %s\n",ha_eq[i].cofname);
    //  for(j=0; j<ha_eq[i].size+1; j++)printf("eq dim %d %s\n",j,ha_eq[i].dimsets[j]);
    //}

    for(i=0; i<neq; i++) {
      /*j3=1;
      j1=0;
      for(j=0; j<ha_eq[i].size; j++)for(j0=0; j0<nset; j0++)if(strcmp(ha_set[ha_eq[i].setid[j]].setname,ha_set[j0].setname)==0) {
            if(j0==allregset) {
              ha_eqreg[i]=j;
              ha_eqtimesbegad[i]=ha_set[j0].begadd;
              j1++;
            }
            ha_eq[i].setid[j]=j0;
            j3*=ha_set[j0].size;
            break;
          }
      if(j1>0) {
        ha_eqint[i]=true;
      }*/
      j3=1;
      if(ha_eqreg[i]>-1)ha_eqtimesbegad[i]=ha_set[ha_eq[i].setid[ha_eqreg[i]]].begadd;
      /*for(j=0; j<ha_eq[i].size; j++){
            j3*=ha_set[ha_eq[i].setid[j]].size;
          }
      ha_eq[i].matsize=j3;*/
      //ha_eqtimematsize[i]=ha_eq[i].matsize;
      ha_eq[i].antidims[ha_eq[i].size-1]=1;
      if(ha_eq[i].size>1) {
        for (j2=ha_eq[i].size-2; j2>-1; j2--) {
          ha_eq[i].antidims[j2]=ha_eq[i].antidims[j2+1]*ha_set[ha_eq[i].setid[j2+1]].size;
        }
      }
    }
    //j1=0;
    //for(i=0; i<neq; i++) j1+=ha_eq[i].matsize;
    //printf("Number of equations %d\n",j1);
    //for(i=0; i<neq; i++)if(!ha_eqint[i])printf("eq %s\n",ha_eq[i].cofname);
    j3=0;
    for (i=0; i<neq; i++) {
      if(ha_eqint[i]) { //for (j=0; j<ha_set[allregset].size; j++)counteq[j]+=(uvadd)ha_eq[i].matsize/ha_set[allregset].size;
        //printf("i %d\n",i);
        for (j=0; j<ha_eq[i].matsize; j++) {
          j0=j;
          for(j1=0; j1<ha_eqreg[i]+1; j1++) {
            j2=j0/ha_eq[i].antidims[j1];
            j0-=j2*ha_eq[i].antidims[j1];
          }
          counteq[ha_setele[ha_eqtimesbegad[i]+j2].setsh[0]]++;
        }
      }
      //j3+=ha_eq[i].matsize;
    }
    //printf("OK!!!\n");
    //for (i=0; i<ha_set[alltimeset].size; i++) if(ha_eqint[i])printf("intraeq %d\n",counteq[i]);
    if(rank==rank_hsl) {
      counteqnoadd[0]=counteq[0];
      j2=counteq[0];
      counteq[0]=0;
      nintraeq=j2;
    }
    for (i=1; i<ndblock; i++) {
      counteqnoadd[i]=counteq[i];
      nintraeq+=counteq[i];
      j3=counteq[i];
      counteq[i]=counteq[i-1]+j2;
      j2=j3;
    }
    if(rank==rank_hsl) {
      counteq[ndblock]=VecSize;
      netcut=VecSize-countvarintra1[ndblock];
      //for (i=0; i<ndblock+1; i++) if(ha_eqint[i])printf("intraeq1 %d total %d\n",counteq[i],nintraeq);
    }
  }
  printf("rank1 %d\n",rank);
  //printf("OK!!!!\n");

  if(alltimeset>=0&&allregset>=0) {
    //if(rank==rank_hsl) {
    //hcge_req(tabfile,commsyntax,ha_eq,neq,ha_set,nset);
    //}
    //for(i=0; i<neq; i++) {
    //  printf("eq name: %s\n",ha_eq[i].cofname);
    //  for(j=0; j<ha_eq[i].size+1; j++)printf("eq dim %d %s\n",j,ha_eq[i].dimsets[j]);
    //}
    switch (nesteddbbd) {
    case 1 :
      for (i=0; i<neq; i++) {
        if(ha_eqtime[i]>-1)ha_eqtimesbegad[i]=ha_set[ha_eq[i].setid[ha_eqtime[i]]].begadd;
        if(ha_eqreg[i]>-1)ha_eqregsbegad[i]=ha_set[ha_eq[i].setid[ha_eqreg[i]]].begadd;
        if(ha_eqint[i]) {
          for (j=0; j<ha_set[ha_eq[i].setid[ha_eqtime[i]]].size; j++)
            //if(ha_eq[i].setid[ha_eqtime[i]]==alltimeset)
            if(ha_eqreg[i]>-1)for(j1=0; j1<ha_set[ha_eq[i].setid[ha_eqreg[i]]].size; j1++)
                counteq[ha_setele[ha_eqtimesbegad[i]+j].setsh[ha_set[ha_eq[i].setid[ha_eqtime[i]]].intsup]*(nreg+1)+ha_setele[ha_eqregsbegad[i]+j1].setsh[ha_set[ha_eq[i].setid[ha_eqreg[i]]].regsup]]+=ha_eq[i].matsize/ha_set[ha_eq[i].setid[ha_eqtime[i]]].size/ha_set[ha_eq[i].setid[ha_eqreg[i]]].size;
            else counteq[ha_setele[ha_eqtimesbegad[i]+j].setsh[ha_set[ha_eq[i].setid[ha_eqtime[i]]].intsup]*(nreg+1)+nreg]+=ha_eq[i].matsize/ha_set[ha_eq[i].setid[ha_eqtime[i]]].size;
        }
      }
      if(rank==rank_hsl) {
        counteqnoadd[0]=counteq[0];
        j2=counteq[0];
        counteq[0]=0;
        nintraeq=j2;
      }
      for (i=1; i<ndblock; i++) {
        counteqnoadd[i]=counteq[i];
        nintraeq+=counteq[i];
        j3=counteq[i];
        counteq[i]=counteq[i-1]+j2;
        j2=j3;
        //printf("counteq %d noadd %d\n",counteq[i],counteqnoadd[i]);
      }
      if(rank==rank_hsl) {
        counteq[ndblock]=VecSize;//Attention!!!!!!!!!!!! Different from countvarintra1. Unchanged for not affecting previous method
        netcut=VecSize-countvarintra1[ndblock];
      }
      //printf("counteq %d\n",counteq[i]);
      //printf("netcut %d nintraeq %d\n",netcut,nintraeq);
      break;

    default :
      for(i=0; i<neq; i++) {
        /*j3=1;
        j1=0;
        j4=0;
        for(j=0; j<ha_eq[i].size; j++)for(j0=0; j0<nset; j0++)if(strcmp(ha_set[ha_eq[i].setid[j]].setname,ha_set[j0].setname)==0) {
              if(j0==alltimeset) {
                ha_eqtime[i]=j;
                ha_eqtimesbegad[i]=ha_set[j0].begadd;
                j1++;
              }
              if(j0==allregset) {
                ha_eqreg[i]=j;
                ha_eqregsbegad[i]=ha_set[j0].begadd;
                j1++;
              }
              ha_eq[i].setid[j]=j0;
              j3*=ha_set[j0].size;
              //break;
            }
        if(j1>1) {
          ha_eqint[i]=true;
        }*/
        //j3=1;
        //j1=0;
        //j4=0;
        if(ha_eqtime[i]>-1)ha_eqtimesbegad[i]=ha_set[ha_eq[i].setid[ha_eqtime[i]]].begadd;
        if(ha_eqreg[i]>-1)ha_eqregsbegad[i]=ha_set[ha_eq[i].setid[ha_eqreg[i]]].begadd;
        /*for(j=0; j<ha_eq[i].size; j++){
              j3*=ha_set[ha_eq[i].setid[j]].size;
            }
        ha_eq[i].matsize=j3;*/
        //ha_eqtimematsize[i]=ha_eq[i].matsize;
//      ha_eq[i].antidims[ha_eq[i].size-1]=1;
//      if(ha_eq[i].size>1) {
//        for (j2=ha_eq[i].size-2; j2>-1; j2--) {
//          ha_eq[i].antidims[j2]=ha_eq[i].antidims[j2+1]*ha_set[ha_eq[i].setid[j2+1]].size;
//        }
//      }
      }
      //j1=0;
      //for(i=0; i<neq; i++) j1+=ha_eq[i].matsize;
      //printf("Number of equations %d\n",j1);
      //j3=0;
      for (i=0; i<neq; i++) {
        if(ha_eqint[i]) {
          for (j=0; j<ha_set[ha_eq[i].setid[ha_eqtime[i]]].size; j++)
            /*for (j=0; j<ha_eq[i].matsize; j++) {
              j0=j;
              for(j1=0; j1<ha_eqtime[i]+1; j1++) {
                j2=j0/ha_eq[i].antidims[j1];
                j0-=j2*ha_eq[i].antidims[j1];
              }
              j0=j;
              for(j1=0; j1<ha_eqreg[i]+1; j1++) {
                j4=j0/ha_eq[i].antidims[j1];
                j0-=j4*ha_eq[i].antidims[j1];
              }
              counteq[ha_setele[ha_eqtimesbegad[i]+j2].setsh*ha_eq[i].dims[ha_eqreg[i]]+ha_setele[ha_eqregsbegad[i]+j4].setsh]++;
            }*/
            if(ha_eq[i].setid[ha_eqtime[i]]==alltimeset)
              //counteq[ha_setele[ha_eqtimesbegad[i]+j].setsh[0]*ha_set[ha_eq[i].setid[ha_eqreg[i]]].size+ha_setele[ha_eqregsbegad[i]+j1].setsh[0]]+=ha_eq[i].matsize/ha_set[ha_eq[i].setid[ha_eqtime[i]]].size/ha_set[ha_eq[i].setid[ha_eqreg[i]]].size;
              for(j1=0; j1<ha_set[ha_eq[i].setid[ha_eqreg[i]]].size; j1++)
                counteq[ha_setele[ha_eqtimesbegad[i]+j].setsh[0]*ha_set[ha_eq[i].setid[ha_eqreg[i]]].size+ha_setele[ha_eqregsbegad[i]+j1].setsh[0]]+=ha_eq[i].matsize/ha_set[ha_eq[i].setid[ha_eqtime[i]]].size/ha_set[ha_eq[i].setid[ha_eqreg[i]]].size;
            else {
              for(j3=1; j3<MAXSUPSET; j3++)if(ha_set[ha_eq[i].setid[ha_eqtime[i]]].subsetid[j3]=alltimeset)break;
              for(j1=0; j1<ha_set[ha_eq[i].setid[ha_eqreg[i]]].size; j1++)
                counteq[ha_setele[ha_eqtimesbegad[i]+j].setsh[j3]*ha_set[ha_eq[i].setid[ha_eqreg[i]]].size+ha_setele[ha_eqregsbegad[i]+j1].setsh[0]]+=ha_eq[i].matsize/ha_set[ha_eq[i].setid[ha_eqtime[i]]].size/ha_set[ha_eq[i].setid[ha_eqreg[i]]].size;
            }
        }
        //j3+=ha_var[i].matsize;
      }
      //for (i=0; i<8*11; i++) printf("intraeq %d\n",counteq[i]);
      if(rank==rank_hsl) {
        counteqnoadd[0]=counteq[0];
        j2=counteq[0];
        counteq[0]=0;
        nintraeq=j2;
      }
      for (i=1; i<ndblock; i++) {
        counteqnoadd[i]=counteq[i];
        nintraeq+=counteq[i];
        j3=counteq[i];
        counteq[i]=counteq[i-1]+j2;
        j2=j3;
      }
      if(rank==rank_hsl) {
        counteq[ndblock]=VecSize;
        netcut=VecSize-countvarintra1[ndblock];
      }
    }
  }
  //for (i=0; i<neq; i++)printf("intraeq1 %d eqtime %d eqreg %d\n",ha_eqint[i],ha_eqtime[i],ha_eqreg[i]);
  if(rank==rank_hsl) {
    //for (i=0; i<ndblock+1; i++)printf("intraeq1 %d total %d\n",counteq[i],nintraeq);
    printf("netcut %ld nintraeq %ld\n",netcut,nintraeq);
  }
  free(ha_eq);
  free(ha_eqtimesbegad);
  free(ha_eqregsbegad);
  //for (i=0; i<nvar; i++)for (j=0; j<ha_var[i].matsize; j++) printf("varele %s exo %d exoindx %ld val %lf\n",ha_var[i].cofname,ha_cgeshock[ha_var[i].begadd+j].ShockId,ha_cgeshock[ha_var[i].begadd+j].ExoIndx,ha_cgeshock[ha_var[i].begadd+j].ShockVal);

  //if(rank==rank_hsl) for (i=0; i<ndblock+1; i++) printf("c %d c1 %d intreq %d\n",counteq[i],countvarintra1[i],nintraeq);

  if(nohsl) {
    VecCreate(PETSC_COMM_WORLD,&vece);
  }
  else {
    VecCreate(PETSC_COMM_SELF,&vece);
  }
  if(nohsl) {
    VecSetType(vece,VECMPI);
  }
  else {
    VecSetType(vece,VECSEQ);
  }
  int localsize=0,nmatint,localbeg,localend;
  int *locals= (int *) calloc (mpisize,sizeof(int));
  if(nesteddbbd==1) {
    nmatint=ntime/mpisize;
    for(i=0; i<mpisize; i++)locals[i]=nmatint;
    for(i=0; i<mpisize; i++)if(i<ntime-mpisize*nmatint)locals[i]++;
    localbeg=0;
    for(i=0; i<mpisize; i++)if(i<rank)localbeg+=locals[i]*(nreg+1);
    localend=0;
    for(i=0; i<mpisize; i++)if(i<rank+1)localend+=locals[i]*(nreg+1);
    if(rank==mpisize-1)localend=ndblock;
    printf("rank %d localbeg %d localend %d\n",rank,localbeg,localend);
    localsize=0;
    for (i=1; i<ndblock+1; i++)if(i>localbeg&&i<=localend)localsize+=counteq[i]-counteq[i-1];
    printf("rank %d localsize %d\n",rank,localsize);
    VecSetSizes(vece,localsize,VecSize);
    //VecSetSizes(vecb,PETSC_DECIDE,VecSize);
    //return 0;
  }
  else {
    VecSetSizes(vece,PETSC_DECIDE,VecSize);
  }
  VecSetOption(vece, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
  free(locals);
  //return 0;
  //ierr = VecDuplicate(x,&b1);
  //CHKERRQ(ierr);
  if(rank==rank_hsl) {
    //VecSet(x,zero);
    //CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(vece,&Istart,&Iend);
    CHKERRQ(ierr);
    //localsize=Iend-Istart;
    ierr = PetscMalloc((Iend-Istart)*sizeof(PetscInt),&dnnz);
    CHKERRQ(ierr);
    ierr = PetscMalloc((Iend-Istart)*sizeof(PetscInt),&onnz);
    CHKERRQ(ierr);
    ierr = PetscMalloc((Iend-Istart)*sizeof(PetscInt),&dnnzB);
    CHKERRQ(ierr);
    ierr = PetscMalloc((Iend-Istart)*sizeof(PetscInt),&onnzB);
    CHKERRQ(ierr);
    for (i=Istart; i<Iend; i++) {
      dnnz[i-Istart]=1;
      onnz[i-Istart]=0;
      dnnzB[i-Istart]=1;
      onnzB[i-Istart]=0;
    }
  }
  //printf("OK");

  printf("rank11 %d Istart %d I end %d\n",rank, Istart,Iend);
  if(rank==rank_hsl) {
    NewMatreadele(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,nexo,ha_cgeshock,ndblock,alltimeset,allregset,ha_eqint,ha_eqadd,ha_eqtime,ha_eqreg,counteq,nintraeq,&sbbd_overrid,Istart,Iend,&dnz,dnnz,&onz,onnz,&dnzB,dnnzB,&onzB,onnzB,nesteddbbd);
  }
  printf("OKla!!!\n");
  if(sbbd_overrid&&!sbbd_overuser) {
    printf(" It looks like you have an intertemporal model, \n please do -enable_time for more accurate results if \n you know what you are doing!\n");
  }
  if(sbbd_overuser) {
    sbbd_overrid=false;
  }
  //for (i=0; i<neq; i++)printf("intraeq1 %d eqtime %d eqreg %d\n",ha_eqint[i],ha_eqtime[i],ha_eqreg[i]);
  //printf("OK!!!\n");
  free(ha_eqint);
  //printf("OK!!!\n");
  free(ha_eqtime);
  //printf("OK!!!\n");
  free(ha_eqreg);
  /*FILE *fout;
  fout = fopen("eq.txt","w");
  for(i=0;i<VecSize;i++){
    fprintf(fout,"%d\n",ha_eqadd[i]);//ha_cgeshock[i].ExoIndx
  }
  fclose(fout);*/
  //printf("nintraeq %d\n",nintraeq);
  //for(i=0;i<ndblock;i++)printf("counteq %d\n",counteqnoadd[i]);
  //printf("OK!!!\n");
  //PetscSynchronizedPrintf(PETSC_COMM_SELF,"Istart: %d Iend %d\n",Istart,Iend);
  //for (i=0; i<10; i++) if(Istart<=i&&i<Iend) PetscPrintf(PETSC_COMM_SELF,"dnnz: %d dnnzB %d i %d\n",dnnz[i-Istart],dnnzB[i-Istart],i);//printf("dnnz %d\n",dnnz[i-Istart]);
  //for (i=0; i<10; i++) if(Istart<=i&&i<Iend) PetscPrintf(PETSC_COMM_SELF,"onnz: %d onnzB %d i %d\n",onnz[i-Istart],onnzB[i-Istart],i);//printf("onnz %d\n",onnz[i-Istart]);
  //PetscPrintf(PETSC_COMM_SELF,"dnz1 %d onz1 %d\n",dnz,onz);

  //ierr = MatCreate(PETSC_COMM_WORLD,&A);
  //CHKERRQ(ierr);
  //ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
  //CHKERRQ(ierr);
  //ierr = MatSetFromOptions(A);
  //CHKERRQ(ierr);
  //ierr = MatMPIAIJSetPreallocation(A,dnz,dnnz,onz,onnz);
  //CHKERRQ(ierr);
  //ierr = MatSeqAIJSetPreallocation(A,dnz+onz,PETSC_NULL);
  //CHKERRQ(ierr);
  //PetscPrintf(PETSC_COMM_SELF,"Ist1!\n");
  dnz=0;
  dnzB=0;
  onz=0;
  onzB=0;
  for (i=Istart; i<Iend; i++) {
//  if(Istart<=i&&i<Iend){
    if (dnnzB[i-Istart]+onnzB[i-Istart]>nexo-1&&dnnzB[i-Istart]>1) {
      dnnzB[i-Istart]--;
    }
    if (dnnz[i-Istart]+onnz[i-Istart]>nvarele-nexo&&dnnz[i-Istart]>1) {
      dnnz[i-Istart]--;
    }
    if (dnnz[i-Istart]>dnz) {
      dnz=dnnz[i-Istart];
    }
    if (onnz[i-Istart]>dnz) {
      onz=dnnz[i-Istart];
    }
    if (dnnzB[i-Istart]>dnzB) {
      dnzB=dnnzB[i-Istart];
    }
    if (onnzB[i-Istart]>onzB) {
      onzB=onnzB[i-Istart];
    }
    //if(i>=5819901&&i<=5819905)printf("i %d dnz %d onz %d\n",i,dnnz[i-Istart],onnz[i-Istart]);
    //if(i==147882)PetscPrintf(PETSC_COMM_SELF,"Ist %d ddn %d onz %d nexo %d dnnzb %d onzb %d\n",i,dnnz[i-Istart],onnz[i-Istart],nexo,dnnzB[i-Istart],onnzB[i-Istart]);
//  }
  }
  printf("OK1!!! rank %d\n",rank);
  //PetscPrintf(PETSC_COMM_SELF,"Ist2!\n");
  //ierr = MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize,dnz,dnnz,onz,onnz,&A);
  //ierr = MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize,0,PETSC_NULL,0,PETSC_NULL,&A);
  //CHKERRQ(ierr);
  /* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
  FILE* tempvar;
  MPI_Fint fcomm;
  fcomm = MPI_Comm_c2f(PETSC_COMM_WORLD);
  PetscInt nz01,*ai,*aj;
  PetscScalar *vals;
  FILE *ofp;
  int sol;
//  struct point indata;
//  indata.m=VecSize;
//  indata.mpisize=mpisize;
//  indata.nblock=ndblock;
//  indata.nsbbdblocks=nsbbdblocks;
//  struct point *ptx= &indata;

  forint indata[5];
  //indata[0]=0;
  indata[1]=VecSize;
  indata[2]=mpisize;
  indata[3]=ndblock;
  indata[4]=nsbbdblocks;
  forint *ptx=NULL;
  ptx = indata;

  ha_cgetype *x1=NULL;//= (ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));
  ha_cgetype *xcf=NULL;
  ha_cgetype *x0=NULL;// (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
  ha_cgetype *b1=NULL;//= (ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));
  //extern void spec51_rank_(int *INSIZE,ha_cgetype *cntl6,int *IRN, int *JCN, ha_cgetype *VA);
  extern void spec48_ssol_(forint *INSIZE,forint *IRN, forint *JCN, ha_cgetype *VA, ha_cgetype *B, ha_cgetype *X);
  extern void spec48_ssol2la_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, ha_cgetype *B, ha_cgetype *X);
  //extern void spec48_msol_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, ha_cgetype *B, ha_cgetype *X, PetscInt *IRNC, PetscInt *JCNC, PetscScalar *VAC,int *IRNB,int *JCNB,PetscScalar *VALUESB,ha_cgetype *VECBIVI,int *bivinzrow0,int *bivinzcol0);//, forint *IRNV, forint *JCNV, ha_cgetype *VAV
  //extern void spec48_esol_(int *INSIZE,int *IRN, ha_cgetype *VA,int *KEEP, ha_cgetype *B, ha_cgetype *X);
  //extern void spec48_rpesol_(int *INSIZE,int *IRN, ha_cgetype *VA,int *KEEP, ha_cgetype *B, ha_cgetype *X,ha_cgetype *cntl,ha_cgetype *rinfo,ha_cgetype *error1,int *icntl,int *info,int *w,int *iw);
  extern void spec48_single_(forint *indata,int *irn, int *jcn,ha_cgetype *b1, ha_cgetype *values,ha_cgetype *x1, int *neleperrow,int *ai1, MPI_Fint *fcomm);
  extern void spec48_nomc66_(forint *indata, int *jcn,ha_cgetype *b1, ha_cgetype *values,ha_cgetype *x1, int *neleperrow, MPI_Fint *fcomm,forint *rowptrin, forint *colptrin);
  //extern void my_spar_add_(ha_cgetype *vecbivi, int *biviindx,int *nz1,ha_cgetype *vecbivi0,int *biviindx0,int *nz0,int *nz2);
  //extern void my_spar_addl_(ha_cgetype *vecbivi, long int *biviindx,int *nz1,ha_cgetype *vecbivi0,long int *biviindx0,int *nz0,int *nz2);
  //extern void my_spar_add1_(ha_cgetype *vecbivi, int *biviindx,int *irn, int *jcn,int *nz1,ha_cgetype *vecbivi0,int *biviindx0,int *nz0,int *nz2,int *ncol);
  //extern void my_spar_add1l_(ha_cgetype *vecbivi, long int *biviindx,int *irn, int *jcn,int *nz1,ha_cgetype *vecbivi0,long int *biviindx0,int *nz0,int *nz2,int *ncol);
  //extern void my_spar_add2_(ha_cgetype *vecbivi, int *biviindx,int *irn, int *jcn,int *nz1,ha_cgetype *vecbivi0,int *biviindx0,int *nz0,int *nz2,int *ncol,ha_cgetype *vecbivi2,int *irn2, int *jcn2,int *j2,ha_cgetype *cntl3);
  //extern void my_spar_comp_(int *biviindx,int *nz1,int *biviindx0,int *nz0,int *nz2);
  //extern void my_spar_compl_(long int *biviindx,int *nz1,long int *biviindx0,int *nz0,int *nz2);
  extern void my_vec_comz_(ha_cgetype vecbivi,int *biviindx,int *col, int *row, int *colsize,int *nz0,int *nz1);
  //extern void prep48_alu_(int *INSIZE,int *IRN,int* JCN,ha_cgetype *VA);
  //extern void prep48_msol_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, PetscInt *IRNC, PetscInt *JCNC, PetscScalar *VAC,int *IRNB,int *JCNB,PetscScalar *VALUESB,ha_cgetype *VECBIVI,long int *bivinzrow0,int *bivinzcol0);//, forint *IRNV, forint *JCNV, ha_cgetype *VAV , ha_cgetype *B
  //extern void spar_mulmin_(ha_cgetype* sol,int* nrow,int* nz,int* irn,int* jcn,ha_cgetype* va,ha_cgetype* res);
  //extern void spar_muladd_(ha_cgetype* sol,int* nrow,int* nz,int* irn,int* jcn,ha_cgetype* va,ha_cgetype* res);
  //extern void spar_vbiviadd_(ha_cgetype* sol,int* bvcol,long int* bvrow,long int* bvsize,int* nrow,int *ncol,int* nz,int* irn,int* jcn,ha_cgetype* va,ha_cgetype* res);

  forint k=0,m=1;
  ha_cgetype temp1,temp2;
  forint tindx1;//,tindx2;
  //debug(&ha_cofvar,&ha_cgeshock,ncofele,nvarele);
  //printf("v %d p %d\n",(*(&ha_cgeshock))[10].ExoIndx,ha_cgeshock[10].ExoIndx);
  //printf("v %lf p %lf\n",(*(&ha_cofvar))[10].varval,ha_cofvar[10].varval);
  printf("rank %d ncof %ld\n",rank,ncof);
  
  if(solmethod==11)Johansen(nohsl,VecSize,A,dnz,dnnz,onz,onnz,B,dnzB,dnnzB,onzB,onnzB,vecb,vece,rank,rank_hsl,mpisize,tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,&ha_cofvar,ncofele+nvarele,ncofele,nvarele,&ha_cgeshock,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,ha_eqadd,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,&xcf);

  if(solmethod==10) { //Johansen

    if(nohsl) {
      MatCreate(PETSC_COMM_WORLD,&A);
    }
    else {
      MatCreate(PETSC_COMM_SELF,&A);
    }
    if(nesteddbbd==1)MatSetSizes(A,localsize,localsize,VecSize,VecSize);
    else MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
    if(nohsl) {
      MatSetType(A,MATMPIAIJ);
    }
    else {
      MatSetType(A,MATSEQAIJ);
    }
    printf("OK1 nohsl %d!!!\n",nohsl);
    if(nohsl) {
      MatMPIAIJSetPreallocation(A,dnz,dnnz,onz,onnz);
    }
    else {
      MatSeqAIJSetPreallocation(A,dnz,dnnz);
    }

    ierr = MatSetOption(A,MAT_SYMMETRIC,PETSC_FALSE);
    CHKERRQ(ierr);
    //ierr = MatSetFromOptions(A);
    //CHKERRQ(ierr);
    //PetscPrintf(PETSC_COMM_SELF,"Ist3!\n");

    //ierr = MatCreate(PETSC_COMM_WORLD,&B);
    //CHKERRQ(ierr);
    //ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,VecSize,nexo);
    //CHKERRQ(ierr);
    //ierr = MatSetFromOptions(B);
    //CHKERRQ(ierr);
    //ierr = MatMPIAIJSetPreallocation(B,dnzB,dnnzB,onzB,onnzB);
    //CHKERRQ(ierr);
    //ierr = MatSeqAIJSetPreallocation(B,dnzB+onzB,PETSC_NULL);
    //CHKERRQ(ierr);

    //ierr = MatGetOwnershipRange(B,&Istart,&Iend);
    //CHKERRQ(ierr);
    //ierr = MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize,dnzB,dnnzB,onzB,onnzB,&B);//nexo
    //ierr = MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,VecSize,nexo,0,PETSC_NULL,0,PETSC_NULL,&B);
    //CHKERRQ(ierr);
    //ierr = MatSetFromOptions(B);
    //CHKERRQ(ierr);
    if(nohsl) {
      MatCreate(PETSC_COMM_WORLD,&B);
    }
    else {
      MatCreate(PETSC_COMM_SELF,&B);
    }
    if(nesteddbbd==1)MatSetSizes(B,localsize,localsize,VecSize,VecSize);
    else MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
    if(nohsl) {
      MatSetType(B,MATMPIAIJ);
    }
    else {
      MatSetType(B,MATSEQAIJ);
    }
    if(nohsl) {
      MatMPIAIJSetPreallocation(B,dnzB,dnnzB,onzB,onnzB);
    }
    else {
      MatSeqAIJSetPreallocation(B,dnzB,dnnzB);
    }
    //ierr = MatGetOwnershipRange(A,&Istart,&Iend);
    //CHKERRQ(ierr);

    //ierr = PetscGetCPUTime(&time1);
    //CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix preparation time %f\n",time1-time0);
    //CHKERRQ(ierr);
    //ierr = PetscGetCPUTime(&time0);
    //CHKERRQ(ierr);
    //printf("rank %d rankhsl %d nohsl %d\n",rank,rank_hsl,nohsl);
    gettimeofday(&endtime, NULL);
    if(rank==0)printf("Matrix preparation time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
    
    if(rank==rank_hsl) {
      //if(solmethod==10||solmethod==0) {
      HaNewMatVal(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,ndblock,alltimeset,allregset,ha_eqadd,counteq,nintraeq,A,B);
      //}
    }
//  if(nohsl) {
//    MPI_Bcast(counteq,(ndblock+1)*sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
//    MPI_Bcast(counteqnoadd,ndblock*sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
//    MPI_Bcast(countvarintra1,(ndblock+1)*sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
//    //MPI_Bcast(ha_eqadd,VecSize*sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
//    //MPI_Bcast(&nexo1,sizeof(uvadd), MPI_BYTE,0, PETSC_COMM_WORLD);
//  }
    //for (i=0; i<ndblock; i++) if(ha_eqint[i])printf("intraeq %d total %d\n",counteq[i],nintraeq);

    //ierr = PetscGetCPUTime(&time1);
    //CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix calculation time %f\n",time1-time0);
    //CHKERRQ(ierr);
    //ierr = PetscGetCPUTime(&time0);
    //CHKERRQ(ierr);
    gettimeofday(&begintime, NULL);
    if(rank==0)printf("Matrix calculation time %f\n",(begintime.tv_sec - endtime.tv_sec)+((double)(begintime.tv_usec - endtime.tv_usec))/ 1000000);
    MPI_Barrier(PETSC_COMM_WORLD);
    if(rank==rank_hsl) {
      strcpy(tempfilenam,temdir);
      strcat(tempfilenam,"_tempvar");
      sprintf(tempchar, "%d",rank);
      strcat(tempfilenam,tempchar);
      strcat(tempfilenam,".bin");
      if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
        printf("Error opening file\n");
        //return 1;
      }
      fwrite(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
      fclose(tempvar);
      free(ha_cofvar);
      ha_cofvar=NULL;
      //ha_cofvar=realloc (ha_cofvar,1*sizeof(ha_cgevar));
    }
    for (count=0; count<nvarele; count++) {
      if (ha_cgeshock[count].ShockId) {
        value = ha_cgeshock[count].ShockVal;
        dnz=ha_cgeshock[count].ExoIndx;
        //printf("dnz %d value %f\n",dnz,value);
        VecSetValues(vece,1,&dnz,&value,INSERT_VALUES);
      }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    ierr = VecAssemblyBegin(vece);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vece);
    CHKERRQ(ierr);
    if(rank==rank_hsl) {
      strcpy(tempfilenam,temdir);
      strcat(tempfilenam,"_tempshock");
      sprintf(tempchar, "%d",rank);
      strcat(tempfilenam,tempchar);
      strcat(tempfilenam,".bin");
      if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
        printf("Error opening file\n");
        //return 1;
      }
      fwrite(ha_cgeshock, sizeof(ha_cgeexovar),nvarele, tempvar);
      fclose(tempvar);
      free(ha_cgeshock);
      ha_cgeshock=NULL;
      //ha_cgeshock=realloc (ha_cgeshock,1*sizeof(ha_cgeexovar));
    }

    //ierr = MatDiagonalSet(A,x,ADD_VALUES);
    //CHKERRQ(ierr);
    MPI_Barrier(PETSC_COMM_WORLD);

    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    printf("OK11\n");
    ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    //ierr = PetscGetCPUTime(&time1);
    //CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix assembly time %f\n",time1-time0);
    //CHKERRQ(ierr);
    //ierr = PetscGetCPUTime(&time0);
    gettimeofday(&endtime, NULL);
    if(rank==0)printf("Matrix assembly time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
    CHKERRQ(ierr);
    PetscViewer viewer;
    //PetscViewerASCIIOpen(PETSC_COMM_WORLD, "A.txt", &viewer);
    //ierr = MatView(A,viewer);
    //PetscViewerDestroy(&viewer);
    //MatView(A,0);
    //neq=ha_cgematdim(tabfile,commsyntax,ha_set,nset);
    //printf("neq %d\n",neq);
    //ha_cgespamat *ha_spamat= (ha_cgespamat *) calloc (nvarele*neq,sizeof(ha_cgespamat));
    //nspamele=ha_cgematrix(tabfile,iodata,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ha_spamat,nvarele);
    //printf("neq %d\n",nspamele);
    //for (i=0;i<nspamele;i++) printf("row %d col %d val %lf\n",ha_spamat[i].row,ha_spamat[i].col,ha_spamat[i].mval);
    //ierr = VecAssemblyBegin(b);
    //CHKERRQ(ierr);
    //ierr = VecAssemblyEnd(b);
    //CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&time1);
    CHKERRQ(ierr);
    //PetscViewerASCIIOpen(PETSC_COMM_WORLD, "B.txt", &viewer);
    //ierr = MatView(B,viewer);
    //CHKERRQ(ierr);
    //MatView(B,0);
    printf("Vector e:\n");
    //VecView(e,0);
    ierr = VecDuplicate(vece,&vecb);
    CHKERRQ(ierr);
    ierr = MatMult(B,vece,vecb);
    CHKERRQ(ierr);
    printf("Vector b:\n");
    ierr = VecAssemblyBegin(vecb);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vecb);
    CHKERRQ(ierr);
    ierr = MatDestroy(&B);
    CHKERRQ(ierr);
    ierr = VecDestroy(&vece);
    CHKERRQ(ierr);
    printf("Vector c:\n");
    //Vec b1,e1;
    //ierr = VecDuplicate(x,&b1);
    //CHKERRQ(ierr);
    //ierr = VecDuplicate(x,&e1);
    //CHKERRQ(ierr);
    //ierr = VecSet(e1,1);
    //CHKERRQ(ierr);
    //ierr = MatMult(A,e1,b1);
    //CHKERRQ(ierr);
    //printf("Vector b1:\n");
    //VecView(b,0);
    //PetscViewer viewer;

    if(matsol>1) {
      //PetscViewerDrawOpen(PETSC_COMM_WORLD,0,"",0,0,700,700,&viewer);
      //ierr = MatView(A,viewer);
      //CHKERRQ(ierr);
      //PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
      //PetscViewerDestroy(&viewer);
      //strcpy(tempchar,"lrzip -D -l -L 1 ");
      //strcat(tempchar,tempfilenam);
      //system(tempchar);
      //remove(tempfilenam);
      //ha_cofvar=realloc (ha_cofvar,sizeof(ha_cgevar));
      //time(&timestr);
      gettimeofday(&begintime, NULL);
      //clock_gettime(CLOCK_REALTIME, &gettime_beg);
      ierr = PetscGetCPUTime(&time1);
      //start_time = gettime_now.tv_nsec;
      int *ha_rows= (int *) calloc (VecSize,sizeof(int));
      int *ha_cols= (int *) calloc (VecSize,sizeof(int));
      int *ha_ndblocks= (int *) calloc (ndblock,sizeof(int));
      //uvadd *orig_exoindx= (uvadd *) calloc (VecSize,sizeof(uvadd));
      if(matsol==2) {
        HaDBBDMatOder(A,VecSize,mpisize,rank,Istart,Iend,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,cntl6);
        x0=realloc (x0,VecSize*sizeof(ha_cgetype));
        //gettimeofday(&endtime, NULL);
        //if(rank==0)printf("Order time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
        //FILE *fout;
        //fout = fopen("eq.txt","w");
        //for(i=0;i<VecSize;i++){
        //fprintf(fout,"%d\n",ha_cols[i]);//ha_cgeshock[i].ExoIndx
        //}
        //fclose(fout);
        //for(i=0;i<ndblock;i++)printf("counteq %d\n",ha_ndblocks[i]);
        //VecView(b,0);
        //VecSet(b,1);
        HaDBBDParSol(A,vecb,x0,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laD,cntl3);//,iter
      }
      if(matsol==3) {
        presol=1;
//        uvadd *counteqs= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
//        uvadd *counteqnoadds= (uvadd *) calloc (ndblock,sizeof(uvadd));
//        uvadd *countvarintra1s= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
        if(presol){
//        memcpy(counteqs,counteq,(ndblock+1)*sizeof(uvadd));
//        memcpy(counteqnoadds,counteqnoadd,(ndblock)*sizeof(uvadd));
//        memcpy(countvarintra1s,countvarintra1,(ndblock+1)*sizeof(uvadd));
        HaNDBBDMatOderPre(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
        //HaDBBDMatOder1(A,VecSize,mpisize,rank,Istart,Iend,ha_cgeshock,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,cntl6);
//        x0=realloc (x0,VecSize*sizeof(ha_cgetype));
        //printf("Order OK!!!\n");
        HaNDBBDParPre(A,vecb,x0,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,nreg,ntime,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
//        printf("Would you like to record ranks?\n");
//        ch = getchar();
        }
        presol=0;
        //memcpy(counteqs,counteq,(ndblock+1)*sizeof(uvadd));
        //memcpy(counteqnoadds,counteqnoadd,(ndblock)*sizeof(uvadd));
        //memcpy(countvarintra1s,countvarintra1,(ndblock+1)*sizeof(uvadd));
        HaNDBBDMatOder(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
        x0=realloc (x0,VecSize*sizeof(ha_cgetype));
        HaNDBBDParSol(A,vecb,x0,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,nreg,ntime,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
//        free(counteqs);
//        free(counteqnoadds);
//        free(countvarintra1s);
      }
      //time(&timeend);
      //ierr = PetscGetCPUTime(&time1);
      //ierr = PetscPrintf(PETSC_COMM_WORLD,"One step solution %f\n",time1-time0);
      //if(rank==0)printf("One step calculation time %f\n",difftime(timeend,timestr));
      gettimeofday(&endtime, NULL);
      //clock_gettime(CLOCK_REALTIME, &gettime_end);
      //rep_time = ((double)(gettime_end.tv_nsec-gettime_beg.tv_nsec))/1000000000.0;
      if(rank==0)printf("One step calculation time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
      //if(rank==0)printf("One step calculation real time %lf\n",rep_time);
      free(ha_rows);
      free(ha_cols);
      free(ha_ndblocks);
      printf("rank %d\n",rank);
      //ierr = MatDestroy(&A);
      //CHKERRQ(ierr);
      //ierr = VecDestroy(&b);
      //CHKERRQ(ierr);
      MPI_Barrier(PETSC_COMM_WORLD);
    }
    else {

      /* This function should be called to be able to use PETSc routines
         from the FORTRAN subroutines needed by this program */
      //PetscInitializeFortran();
      //time(&timestr);
      gettimeofday(&begintime, NULL);
      ierr = PetscGetCPUTime(&time0);
      MPI_Fint fcomm;
      fcomm = MPI_Comm_c2f(PETSC_COMM_WORLD);
      Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)A->data;
      FILE *ofp;

      if(rank==rank_hsl) {
        ai= aa->i;
        aj= aa->j;
        vals=aa->a;
        nz01=aa->nz;
        //for(i=0; i<10; i++) printf("x %f irn %d jcn %d\n",vals[i],ai[i],aj[i]);
        count=0;
        for(i=0; i<nz01; i++) if(vals[i]!=0) {
            count++;
          }
        printf("count %d nz %d\n",count,nz01);
      }
      indata[0]=count;//.nz
      //if(alltimeset>=0||allregset>=0)la=1;
      //MPI_Bcast(&count,1, MPI_LONG, 0, PETSC_COMM_WORLD);
      //MPI_Bcast(&nz,1, MPI_LONG, 0, PETSC_COMM_WORLD);

      if(matsol==1) {
        int *irn=(int *) calloc (count,sizeof(int));
        int *irn1=(int *) calloc (nz01,sizeof(int));
        int *jcn=(int *) calloc (count,sizeof(int));
        ha_cgetype *values= (ha_cgetype *) calloc (count,sizeof(ha_cgetype));
        if(rank==rank_hsl) {
          for(i=0; i<VecSize-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
              irn1[j]=i+1;
            }
          for(j=ai[VecSize-1]; j<nz01; j++) {
            irn1[j]=VecSize;
          }
          //printf("count %d\n",count);
          //for (i=108494341;i<VecSize;i++)for(j=ai[i]; j<ai[i+1]; j++)printf("Errors here i %ld ai %ld j %ld val %lf\n",i,ai[i],j,vals[j]);
          j=0;
          for(i=0; i<nz01; i++) if(vals[i]!=0) {
              irn[j]=irn1[i];
              jcn[j]=aj[i]+1;
              values[j]=vals[i];
              j++;
            }
        }
        free(irn1);
        ierr = MatDestroy(&A);
        CHKERRQ(ierr);

        //ha_cgetype *x2s= (ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));
        //ha_cgetype *x3s= (ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));
        //ha_cgetype *stepsize= (ha_cgetype *) calloc (nexo,sizeof(ha_cgetype));
        b1=realloc (b1,VecSize*sizeof(ha_cgetype));
        if(rank==rank_hsl) {
          VecGetArray(vecb,&vals);
          for(i=0; i<VecSize; i++) {
            b1[i]=vals[i];
          }
        }
        ierr = VecDestroy(&vecb);
        CHKERRQ(ierr);

        int *neleperrow= (int *) calloc (VecSize,sizeof(int));
        int *ai1= (int *) calloc (VecSize,sizeof(int));
        if(rank==rank_hsl) {
          j=1;
          for(i=1; i<count; i++) {
            if(irn[i]-irn[i-1]>0) {
              neleperrow[k]=j;
              ai1[k]=m;
              j=1;
              m=i+1;
              //if(k+1<irn[i-1]){printf("irn %ld i111 %ld k %ld\n",irn[i-1],i,k);return 0;}
              //printf("k %d ai1 %d nele %d count %d\n",k,ai1[k],neleperrow[k],count);
              k++;
            }
            else {
              j++;
            }
            //if(i>25043001-100&&i<25043001+100)printf("irn %ld k %ld j %ld\n",irn[i],k,j);
          }
          neleperrow[k]=j;
          //printf("vec %ld k %ld j %ld\n",VecSize,k,j);
          ai1[k]=ai1[k-1]+neleperrow[k-1];
          //printf("k %d ai1 %d nele %d count %d\n",k,ai1[k],neleperrow[k],count);
          //for(i=0;i<VecSize;i++)printf("nele %d ai1 %d row %d row1 %d\n",neleperrow[i],ai1[i],i,irn[ai1[i]]);
          //for(i=count-100; i<count; i++) printf("x %f irn %d jrn %d count %d nz %d\n",values[i],irn[i],jcn[i],count,nz);
          //VecGetArray(b2,&vals);
        }
        //for (i=0;i<VecSize;i++)if(neleperrow[i]==0)printf("Errors here %ld\n",i);
        //return 0;
        //free(irn1);
        //ptx= &indata;
        //ierr = MatDestroy(&G);
        //CHKERRQ(ierr);
        //ierr = VecDestroy(&b2);
        //CHKERRQ(ierr);
        //PetscViewerASCIIOpen(PETSC_COMM_WORLD, "A1.txt", &viewer);
        //ierr = MatView(A,viewer);
        //CHKERRQ(ierr);
        //PetscViewerDrawOpen(PETSC_COMM_WORLD,0,"",0,0,700,700,&viewer);
        //ierr = MatView(A,viewer);
        //CHKERRQ(ierr);
        //PetscViewerDestroy(&viewer);
        //ierr = MatDestroy(&A);
        //CHKERRQ(ierr);
        //printf("myid %d time %d reg %d\n",rank,alltimeset,allregset);
        //if(alltimeset>=0||allregset>=0) {
        //printf("add %d int %d fint %d\n",sizeof(uvadd),sizeof(int),sizeof(MPI_Fint));
        x0=realloc (x0,VecSize*sizeof(ha_cgetype));
        if(mc66!=0)spec48_single_(ptx,irn,jcn,b1,values,x0,neleperrow,ai1,&fcomm);
        free(irn);
        if(mc66==0)spec48_nomc66_(ptx,jcn,b1,values,x0,neleperrow,&fcomm,counteq,countvarintra1);
        free(jcn);
        free(values);
        free(neleperrow);
        free(ai1);
        free(b1);
        b1=NULL;
        //b1=realloc (b1,sizeof(ha_cgetype));
      }
      else {
        x0=realloc (x0,VecSize*sizeof(ha_cgetype));
        //b1=realloc (b1,VecSize*sizeof(ha_cgetype));
        lasize=ceil((laA/100)*count);
        int *irn=(int *) calloc (lasize,sizeof(int));
        int *irn1=(int *) calloc (nz01,sizeof(int));
        int *jcn=(int *) calloc (lasize,sizeof(int));
        ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
        //forint *neleperrow= (forint *) calloc (VecSize,sizeof(forint));
        //forint *ai1= (forint *) calloc (VecSize,sizeof(forint));

        if(rank==rank_hsl) {
          for(i=0; i<VecSize-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
              irn1[j]=i+1;
            }
          for(j=ai[VecSize-1]; j<nz01; j++) {
            irn1[j]=VecSize;
          }
          //printf("count %d\n",count);
          j=0;
          for(i=0; i<nz01; i++) if(vals[i]!=0) {
              irn[j]=irn1[i];
              jcn[j]=aj[i]+1;
              values[j]=vals[i];
              j++;
            }
          //for(i=0; i<count; i++) printf("irn %d jcn %d count %d nz %d x %f\n",irn[i],jcn[i],count,i,values[i]);
          /*j=1;
          for(i=1; i<count; i++) {
            if(irn[i]-irn[i-1]>0) {
              neleperrow[k]=j;
              ai1[k]=m;
              j=1;
              m=i+1;
              //printf("k %d ai1 %d nele %d count %d\n",k,ai1[k],neleperrow[k],count);
              k++;
            } else {
              j++;
            }
          }
          neleperrow[k]=j;
          ai1[k]=ai1[k-1]+neleperrow[k-1];*/
          //printf("k %d ai1 %d nele %d count %d\n",k,ai1[k],neleperrow[k],count);
          //for(i=0;i<VecSize;i++)printf("nele %d ai1 %d row %d row1 %d\n",neleperrow[i],ai1[i],i,irn[ai1[i]]);
          //for(i=count-100; i<count; i++) printf("x %f irn %d jrn %d count %d nz %d\n",values[i],irn[i],jcn[i],count,nz);
          //VecGetArray(b2,&vals);
          VecGetArray(vecb,&vals);
          //for(i=0; i<VecSize; i++) {
          //b1[i]=vals[i];
          //}
          //for(i=0; i<VecSize; i++)printf("b1 %f\n",b1[i]);
        }
        //return 0;
        free(irn1);
        PetscViewerDrawOpen(PETSC_COMM_WORLD,0,"",0,0,500,500,&viewer);
        ierr = MatView(A,viewer);
        //CHKERRQ(ierr);
        PetscViewerDestroy(&viewer);
        //PetscViewerASCIIOpen(PETSC_COMM_WORLD, "A1.txt", &viewer);
        //ierr = MatView(A,viewer);
        //CHKERRQ(ierr);
        //CHKERRQ(ierr);
        //PetscViewerDestroy(&viewer);
        ierr = MatDestroy(&A);
        CHKERRQ(ierr);
        //free(ha_cofvar);
        //ierr = VecDestroy(&b);
        //CHKERRQ(ierr);
        //MPI_Bcast(irn, count, MPI_LONG, 0, PETSC_COMM_WORLD);
        //MPI_Bcast(jcn, count, MPI_LONG, 0, PETSC_COMM_WORLD);
        //MPI_Bcast(b1, VecSize, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
        //MPI_Bcast(values, count, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
        //MPI_Bcast(neleperrow, VecSize, MPI_LONG, 0, PETSC_COMM_WORLD);
        //MPI_Bcast(ai1, VecSize, MPI_LONG, 0, PETSC_COMM_WORLD);

        int *insize=(int *) calloc (4,sizeof(int));
        insize[0]=VecSize;
        insize[1]=VecSize;
        insize[2]=count;
        insize[3]=laA;
        //spec48_ssol_(insize,irn,jcn,values,b1,x0);
        if(rank==rank_hsl)spec48_ssol2la_(insize,irn,jcn,values,vals,x0);
        free(insize);
        free(irn);
        free(jcn);
        free(values);
        //free(neleperrow);
        //free(ai1);
        free(b1);
        b1=NULL;
        //b1=realloc (b1,sizeof(ha_cgetype));
        ierr = VecDestroy(&vecb);
        CHKERRQ(ierr);
      }
      //time(&timeend);
      ierr = PetscGetCPUTime(&time1);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"One step solution %f\n",time1-time0);
      //if(rank==0)printf("One step calculation time %f\n",difftime(timeend,timestr));
      gettimeofday(&endtime, NULL);
      if(rank==0)printf("One step calculation time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
    }
    if(rank==rank_hsl) {
      if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
        printf("Error opening file\n");
      }
      ha_cgeshock=realloc (ha_cgeshock,(nvarele)*sizeof(ha_cgeexovar));
      freadresult=fread(ha_cgeshock, sizeof(ha_cgeexovar),nvarele, tempvar);
      fclose(tempvar);
      remove(tempfilenam);

      strcpy(tempfilenam,temdir);
      strcat(tempfilenam,"_tempvar");
      sprintf(tempchar, "%d",rank);
      strcat(tempfilenam,tempchar);
      strcat(tempfilenam,".bin");
      if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
        printf("Error opening file\n");
      }
      ha_cofvar=realloc (ha_cofvar,(ncofele+nvarele)*sizeof(ha_cgevar));
      freadresult=fread(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
      fclose(tempvar);
      remove(tempfilenam);
    }
    //for(i=0; i<VecSize; i++)printf("x %f\n",x1[i]);
    //ha_cgetype *xs= (ha_cgetype *) calloc (nvarele+ncofele,sizeof(ha_cgetype));
    //ha_cgetype *xc= (ha_cgetype *) calloc (nvarele,sizeof(ha_cgetype));
    xcf=realloc (xcf,nvarele*sizeof(ha_cgetype));
    printf("Hello world1!\n");
    if(rank==rank_hsl) {
      //ha_cgetype *varchange= (ha_cgetype *) calloc (nvarele,sizeof(ha_cgetype));
      ha_cofvar1=ha_cofvar+ncofele;
      for(i=0; i<nvar; i++) {
        if(ha_var[i].change_real) {
          for(j=ha_var[i].begadd; j<ha_var[i].matsize+ha_var[i].begadd; j++) {
            //tindx2=ncofele+j;
            if(ha_cgeshock[j].ShockId) {
              ha_cofvar1[j].var0=ha_cofvar1[j].varval;
              ha_cofvar1[j].varval+=ha_cgeshock[j].ShockVal;
              xcf[j]=ha_cgeshock[j].ShockVal;//varchange[j]
              ha_cofvar1[j].csolpupd=ha_cgeshock[j].ShockVal;
            }
            else {
              ha_cofvar1[j].var0=ha_cofvar1[j].varval;
              ha_cofvar1[j].varval+=x0[ha_cgeshock[j].ExoIndx];
              xcf[j]=x0[ha_cgeshock[j].ExoIndx];//varchange[j]
              ha_cofvar1[j].csolpupd=x0[ha_cgeshock[j].ExoIndx];
            }
          }
        }
        else {
          for(j=ha_var[i].begadd; j<ha_var[i].matsize+ha_var[i].begadd; j++) {
            //tindx2=ncofele+j;
            if(ha_cgeshock[j].ShockId) {
              ha_cofvar1[j].var0=ha_cofvar1[j].varval;
              ha_cofvar1[j].varval+=ha_cgeshock[j].ShockVal*ha_cofvar1[j].var0/100;
              xcf[j]=ha_cgeshock[j].ShockVal;//varchange[j]
              ha_cofvar1[j].csolpupd=ha_cgeshock[j].ShockVal;
            }
            else {
              ha_cofvar1[j].var0=ha_cofvar1[j].varval;
              xcf[j]=x0[ha_cgeshock[j].ExoIndx];//varchange[j]
              ha_cofvar1[j].varval+=x0[ha_cgeshock[j].ExoIndx]/100*ha_cofvar1[j].varval;
              ha_cofvar1[j].csolpupd=x0[ha_cgeshock[j].ExoIndx];
            }
          }
        }
      }
      printf("Rank %d Hello world1a!\n",rank);
      hnew_update(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele);
      strcpy(commsyntax,"formula");
      IsIni=false;
      printf("Rank %d Hello world1b!\n",rank);
      hnew_calcff(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,IsIni);
      printf("Rank %d Hello world1c!\n",rank);

      /*if(sol==0)for(i=0; i<ncofele+nvarele; i++) {
          xs[i]=ha_cofvar[i].varval;
        }*/
      //if(sol==0)
//      for(i=0; i<nvarele; i++) {
//        xcf[i]=varchange[i];
//      }
//      free(varchange);
    }
    //printf("Rank %d Hello world2!\n",rank);
  }
  
  int stepcount;
  int nsteps=3;
  ha_cgetype vpercents=1.0,perprecis=0;
  //ierr = VecDuplicate(b,&e);
  //CHKERRQ(ierr);
  //uvadd *ha_eqadd1= (uvadd *) calloc (1,sizeof(uvadd));//recycle ha_cgeset
  FILE* solution;
  int maxsol=3;
    if(solmethod==0)ModMidPoint(nohsl,VecSize,&A,dnz,dnnz,onz,onnz,&B,dnzB,dnnzB,onzB,onnzB,&vecb,&vece,rank,rank_hsl,mpisize,tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,&ha_cofvar,ncofele+nvarele,ncofele,nvarele,&ha_cgeshock,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,ha_eqadd,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,subints,fcomm,&xcf,0);

    if(solmethod==1) { //Modified midpoint Pearson 1991
              uvadd *counteqs= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
              uvadd *counteqnoadds= (uvadd *) calloc (ndblock,sizeof(uvadd));
              uvadd *countvarintra1s= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
              memcpy(counteqs,counteq,(ndblock+1)*sizeof(uvadd));
              memcpy(counteqnoadds,counteqnoadd,(ndblock)*sizeof(uvadd));
              memcpy(countvarintra1s,countvarintra1,(ndblock+1)*sizeof(uvadd));
    //ierr = VecDuplicate(b,&e);
    //CHKERRQ(ierr);
    //time(&timemulti);
    gettimeofday(&begintime, NULL);
    ha_cgetype *xc0=NULL;// (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
    //ha_cgetype *xc012= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
    //ha_cgetype *xc024= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
    ha_cgetype *xc12= NULL;//(ha_cgetype *) calloc (1,sizeof(ha_cgetype));
    ha_cgetype *xc24= NULL;//(ha_cgetype *) calloc (1,sizeof(ha_cgetype));
    int *xc124= NULL;//(int *) calloc (1,sizeof(int));
    ha_cgetype *clag1= NULL;//(ha_cgetype *) calloc (nvarele,sizeof(ha_cgetype));
    ha_cgetype *varchange= NULL;//(ha_cgetype *) calloc (nvarele,sizeof(ha_cgetype));
    for(subindx=0; subindx<subints; subindx++) {
      for(sol=0; sol<maxsol; sol++) {
//         if(sol<2) {
//           nsteps=2*(sol+1);
//         }
//         else {
//           nsteps=8;
//         }
        if(sol==0)nsteps=step1;
        if(sol==1) nsteps=(int)step1*kindx1;
        if(sol==2) nsteps=(int)step1*kindx2;
        vpercents=(ha_cgetype)100/nsteps;
        //nsteps=5;
        //nsteps=1;
        for(stepcount=0; stepcount<nsteps; stepcount++) {
          //for (i=0; i<(nvarele+ncofele); i++) printf("%d cofvar %s val %lf\n",i-ncofele,ha_cofvar[i].varname,ha_cofvar[i].varval);
          printf("subint %d sol %d stepcount %d nsteps %d\n",subindx,sol,stepcount,nsteps);
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time0);
          CHKERRQ(ierr);
          if(stepcount==0) {
            //if(subindx>0)for
            //if(!(sol==0&&subindx==0)) {//if(sol!=0) {
            MPI_Barrier(PETSC_COMM_WORLD);
            if(!(subindx==0&&sol==0&&stepcount==0)) {
              if(nohsl) {
                VecCreate(PETSC_COMM_WORLD,&vece);
              }
              else {
                VecCreate(PETSC_COMM_SELF,&vece);
              }
              if(nohsl) {
                VecSetType(vece,VECMPI);
              }
              else {
                VecSetType(vece,VECSEQ);
              }
              if(nesteddbbd==1)VecSetSizes(vece,localsize,VecSize);
              else VecSetSizes(vece,PETSC_DECIDE,VecSize);
              VecSetOption(vece, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
            }
            if(sol==0)for(i=0; i<ncofele; i++) {
                ha_cofvar[i].var0=ha_cofvar[i].varval;
              }
            else for(i=0; i<ncofele; i++) {
                ha_cofvar[i].varval=ha_cofvar[i].var0;
              }
            ha_cofvar1=ha_cofvar+ncofele;
            for(i=0; i<nvar; i++) {
              if(ha_var[i].change_real) {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  //tindx1=ha_var[i].begadd+j;
                  //tindx2=tindx1+ncofele;
                  if(ha_cgeshock[tindx1].ShockId) {
                    if(sol==0) {
                      ha_cofvar1[tindx1].var0=ha_cofvar1[tindx1].varval;
                    }
                    else {
                      ha_cofvar1[tindx1].varval=ha_cofvar1[tindx1].var0;
                    }
                    ha_cofvar1[tindx1].csolpupd=ha_cgeshock[tindx1].ShockVal/nsteps;
                    //if(ha_cofvar[ncofele+ha_var[i].begadd+j].csolpupd!=0) printf("shock %lf val %lf\n",ha_cofvar[ncofele+ha_var[i].begadd+j].csolpupd,ha_cgeshock[ha_var[i].begadd+j].ShockVal);
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,ha_cofvar1[tindx1].csolpupd,INSERT_VALUES);
                  }
                  else {
                    if(sol==0) {
                      ha_cofvar1[tindx1].var0=ha_cofvar1[tindx1].varval;
                    }
                    else {
                      ha_cofvar1[tindx1].varval=ha_cofvar1[tindx1].var0;
                    }
                  }
                  //ha_cofvar[ncofele+ha_var[i].begadd+j].varval=ha_cofvar[ncofele+ha_var[i].begadd+j].varval+x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx]/nsteps;
                  //ha_cofvar[ncofele+ha_var[i].begadd+j].varchange=x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx]/nsteps;
                }
              }
              else {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  //tindx1=ha_var[i].begadd+j;
                  //tindx2=tindx1+ncofele;
                  if(ha_cgeshock[tindx1].ShockId) {
                    if(sol==0) {
                      ha_cofvar1[tindx1].var0=ha_cofvar1[tindx1].varval;
                    }
                    else {
                      ha_cofvar1[tindx1].varval=ha_cofvar1[tindx1].var0;
                    }
                    temp2=ha_cgeshock[tindx1].ShockVal;//subints;
                    ha_cofvar1[tindx1].csolpupd=(100+(subindx+1)*temp2)/(100+subindx*temp2)-1;//ha_cgeshock[ha_var[i].begadd+j].ShockVal/nsteps;//(exp(log(1+ha_cgeshock[ha_var[i].begadd+j].ShockVal/100)/nsteps)-1)*100;
                    ha_cofvar1[tindx1].csolpupd*=vpercents;//nsteps*100;
                    //if(ha_cofvar[ncofele+tindx1].csolpupd>0)printf("shock %lf val %lf v %lf\n",ha_cofvar[ncofele+ha_var[i].begadd+j].csolpupd,ha_cgeshock[ha_var[i].begadd+j].ShockVal,vpercents);
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,ha_cofvar1[tindx1].csolpupd,INSERT_VALUES);
                  }
                  else {
                    if(sol==0) {
                      ha_cofvar1[tindx1].var0=ha_cofvar1[tindx1].varval;
                    }
                    else {
                      ha_cofvar1[tindx1].varval=ha_cofvar1[tindx1].var0;
                    }
                  }
                  //ha_cofvar[ncofele+ha_var[i].begadd+j].varval*=1+x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx]/nsteps/100;
                  //ha_cofvar[ncofele+ha_var[i].begadd+j].varchange=x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx]/nsteps;
                }
              }
            }
            MPI_Barrier(PETSC_COMM_WORLD);
            ierr = VecAssemblyBegin(vece);
            CHKERRQ(ierr);
            ierr = VecAssemblyEnd(vece);
            CHKERRQ(ierr);
          }
          //printf("OK!!!\n");
          //for (i=0; i<10; i++) printf("%d cofvar %s val %lf change %lf\n",i-ncofele,ha_cofvar[i].varname,ha_cofvar[i].varval,ha_cofvar[i].varchange);
          //for (i=0; i<(nvarele+ncofele); i++) printf("%d cofvar1 %s val %lf change %lf\n",i-ncofele,ha_cofvar[i].varname,ha_cofvar[i].varval,ha_cofvar[i].varchange);
          if(rank==rank_hsl) {
//            if(sizeof(xc0)>sizeof(ha_cgetype)) {
//              strcpy(tempfilenam,temdir);
//              strcat(tempfilenam,"_tempxcO");
//              sprintf(tempchar, "%d",rank);
//              strcat(tempfilenam,tempchar);
//              strcat(tempfilenam,".bin");
//              if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
//                printf("Error opening file\n");
//                //return 1;
//              }
//              fwrite(xc0, sizeof(ha_cgetype),nvarele, tempvar);
//              fclose(tempvar);
//              xc0=realloc (xc0,2*sizeof(ha_cgetype));
//            }
//            strcpy(tempfilenam,temdir);
//            strcat(tempfilenam,"_tempxcf");
//            sprintf(tempchar, "%d",rank);
//            strcat(tempfilenam,tempchar);
//            strcat(tempfilenam,".bin");
//            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
//              printf("Error opening file\n");
//              //return 1;
//            }
//            fwrite(xcf, sizeof(ha_cgetype),nvarele, tempvar);
//            fclose(tempvar);
//            xcf=realloc (xcf,1*sizeof(ha_cgetype));

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempclag1");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
              //return 1;
            }
            fwrite(clag1, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            free(clag1);
            clag1=NULL;
            //clag1=realloc (clag1,1*sizeof(ha_cgetype));

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempvarchange");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
              //return 1;
            }
            fwrite(varchange, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            free(varchange);
            varchange=NULL;
            //varchange=realloc (varchange,1*sizeof(ha_cgetype));
          }
          MPI_Barrier(PETSC_COMM_WORLD);

          strcpy(commsyntax,"equation");
          //ierr = MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize,dnz,dnnz,onz,onnz,&A);
          //CHKERRQ(ierr);
          //ierr = MatSetFromOptions(A);
          //CHKERRQ(ierr);
          //ierr = MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize,dnzB,dnnzB,onzB,onnzB,&B);
          //CHKERRQ(ierr);
          //ierr = MatSetFromOptions(B);
          //CHKERRQ(ierr);

          if(nohsl) {
            MatCreate(PETSC_COMM_WORLD,&A);
          }
          else {
            MatCreate(PETSC_COMM_SELF,&A);
          }
          if(nesteddbbd==1)MatSetSizes(A,localsize,localsize,VecSize,VecSize);
          else MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
          if(nohsl) {
            MatSetType(A,MATMPIAIJ);
          }
          else {
            MatSetType(A,MATSEQAIJ);
          }
          if(nohsl) {
            MatMPIAIJSetPreallocation(A,dnz,dnnz,onz,onnz);
          }
          else {
            MatSeqAIJSetPreallocation(A,dnz,dnnz);
          }

          if(nohsl) {
            MatCreate(PETSC_COMM_WORLD,&B);
          }
          else {
            MatCreate(PETSC_COMM_SELF,&B);
          }
          if(nesteddbbd==1)MatSetSizes(B,localsize,localsize,VecSize,VecSize);
          else MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
          if(nohsl) {
            MatSetType(B,MATMPIAIJ);
          }
          else {
            MatSetType(B,MATSEQAIJ);
          }
          if(nohsl) {
            MatMPIAIJSetPreallocation(B,dnzB,dnnzB,onzB,onnzB);
          }
          else {
            MatSeqAIJSetPreallocation(B,dnzB,dnnzB);
          }

          if(rank==rank_hsl) {
            HaNewMatVal(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,ndblock,alltimeset,allregset,ha_eqadd,counteq,nintraeq,A,B);
          }
          //printf("OK!!!\n");
          if(rank==rank_hsl) {
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempvar");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
              //return 1;
            }
            fwrite(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
            fclose(tempvar);
            free(ha_cofvar);
            ha_cofvar=NULL;
            //ha_cofvar=realloc (ha_cofvar,1*sizeof(ha_cgevar));

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempshock");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
              //return 1;
            }
            fwrite(ha_cgeshock, sizeof(ha_cgeexovar),nvarele, tempvar);
            fclose(tempvar);
            free(ha_cgeshock);
            ha_cgeshock=NULL;
            //ha_cgeshock=realloc (ha_cgeshock,1*sizeof(ha_cgeexovar));
          }

          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
          CHKERRQ(ierr);
          ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
          CHKERRQ(ierr);
          ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
          CHKERRQ(ierr);
          ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
          CHKERRQ(ierr);
          //printf("here!!!!!!!!!!!!!!!!\n");
          //if(!(subindx==0&&sol==0&&stepcount==0)) {
          ierr = VecDuplicate(vece,&vecb);
          CHKERRQ(ierr);
          //}
          if(rank==rank_hsl) {
            ierr = MatMult(B,vece,vecb);
            CHKERRQ(ierr);
          }
          ierr = VecDestroy(&vece);
          CHKERRQ(ierr);
          //printf("here1!!!!!!!!!!!!!!!!\n");
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = VecAssemblyBegin(vecb);
          CHKERRQ(ierr);
          ierr = VecAssemblyEnd(vecb);
          CHKERRQ(ierr);
          //printf("here!!!!!!!!!!!!!!!!\n");
          /*PetscViewerBinaryOpen(PETSC_COMM_WORLD, "b.out",FILE_MODE_WRITE,&viewer);
          ierr = VecView(b,viewer);
          CHKERRQ(ierr);
          PetscViewerDestroy(&viewer);
          //PetscViewerASCIIOpen(PETSC_COMM_WORLD, "b1.txt", &viewer);
          //ierr = VecView(b,viewer);
          //CHKERRQ(ierr);
          //PetscViewerDestroy(&viewer);
          VecCreateSeq(PETSC_COMM_SELF,VecSize,&b2);
          if(rank==0) {
            PetscViewerBinaryOpen(PETSC_COMM_SELF,"b.out",FILE_MODE_READ,&viewer);
            VecLoad(b2,viewer);
            PetscViewerDestroy(&viewer);
          }
          PetscViewerBinaryOpen(PETSC_COMM_WORLD, "A1.out",FILE_MODE_WRITE, &viewer);
          ierr = MatView(A,viewer);
          CHKERRQ(ierr);
          PetscViewerDestroy(&viewer);*/
          //Mat Atrans;
          //MatTranspose(A, MAT_INITIAL_MATRIX,&Atrans);
          //PetscViewerASCIIOpen(PETSC_COMM_WORLD, "A1.txt", &viewer);
          //ierr = MatView(A,viewer);
          //CHKERRQ(ierr);
          //PetscViewerDestroy(&viewer);
          /*ierr = MatDestroy(&A);
          CHKERRQ(ierr);*/
          ierr = MatDestroy(&B);
          CHKERRQ(ierr);
          //MatCreate(PETSC_COMM_SELF,&G);
          //MatSetType(G,MATSEQAIJ);
          /*if(rank==rank_hsl) {
            if ( (tempvar = fopen("_tempvar.bin", "wb")) == NULL ) {
              printf("Error opening file\n");
              //return 1;
            }
            fwrite(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
            fclose(tempvar);
            ha_cofvar=realloc (ha_cofvar,sizeof(ha_cgevar));
          }*/
          if(matsol>=2) {
            //strcpy(tempchar,"lrzip -D -l -L 1 ");
            //strcat(tempchar,tempfilenam);
            //system(tempchar);
            //remove(tempfilenam);
            //ha_cofvar=realloc (ha_cofvar,sizeof(ha_cgevar));
            int *ha_rows= (int *) calloc (VecSize,sizeof(int));
            int *ha_cols= (int *) calloc (VecSize,sizeof(int));
            int *ha_ndblocks= (int *) calloc (ndblock,sizeof(int));
            //ha_eqadd1=realloc(ha_eqadd1,VecSize*sizeof(uvadd));
            //memcpy (ha_eqadd1,ha_eqadd,VecSize*sizeof(uvadd));
            //uvadd *orig_exoindx= (uvadd *) calloc (VecSize,sizeof(uvadd));
            //timestr=clock();//time(&timestr);
            gettimeofday(&gettime_now, NULL);

            if(matsol==2) {
              HaDBBDMatOder(A,VecSize,mpisize,rank,Istart,Iend,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,cntl6);
              x1=realloc (x1,VecSize*sizeof(ha_cgetype));
              //for(i=0;i<VecSize;i++)printf("ha_row %d\n",ha_cols[i]);
              //VecView(b,0);
              //VecSet(b,1);
              HaDBBDParSol(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laD,cntl3);//,iter
            }

            if(matsol==3) {
              presol=1;
//              uvadd *counteqs= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
//              uvadd *counteqnoadds= (uvadd *) calloc (ndblock,sizeof(uvadd));
//              uvadd *countvarintra1s= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
//              if(presol){
              memcpy(counteq,counteqs,(ndblock+1)*sizeof(uvadd));
              memcpy(counteqnoadd,counteqnoadds,(ndblock)*sizeof(uvadd));
              memcpy(countvarintra1,countvarintra1s,(ndblock+1)*sizeof(uvadd));
              HaNDBBDMatOderPre(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
              HaNDBBDParPre(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,nreg,ntime,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
//        printf("Would you like to record ranks?\n");
//        ch = getchar();
//              }
              presol=0;
//              memcpy(counteqs,counteq,(ndblock+1)*sizeof(uvadd));
//              memcpy(counteqnoadds,counteqnoadd,(ndblock)*sizeof(uvadd));
//              memcpy(countvarintra1s,countvarintra1,(ndblock+1)*sizeof(uvadd));
              HaNDBBDMatOder(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
              x1=realloc (x1,VecSize*sizeof(ha_cgetype));
              HaNDBBDParSol(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,nreg,ntime,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
//              free(counteqs);
//              free(counteqnoadds);
//              free(countvarintra1s);
            }

            //timeend=clock();//time(&timeend);
            gettimeofday(&endtime, NULL);
            //memcpy (ha_eqadd,ha_eqadd1,VecSize*sizeof(uvadd));
            //ha_eqadd1=realloc(ha_eqadd1,1*sizeof(uvadd));
            MPI_Barrier(PETSC_COMM_WORLD);
            ierr = PetscGetCPUTime(&time1);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"One step solution %f\n",time1-time0);
            //if(rank==0)printf("One step calculation time %f\n",difftime(timeend,timestr)/CLOCKS_PER_SEC);
            if(rank==0)printf("One step calculation time %f\n",(endtime.tv_sec - gettime_now.tv_sec)+((double)(endtime.tv_usec - gettime_now.tv_usec))/ 1000000);
            free(ha_rows);
            free(ha_cols);
            free(ha_ndblocks);
            printf("rank %d\n",rank);
            MPI_Barrier(PETSC_COMM_WORLD);
          }
          else {
            if(matsol==1) {
              if(rank==rank_hsl) {
                //PetscViewerBinaryOpen(PETSC_COMM_SELF,"A1.out",FILE_MODE_READ,&viewer);
                //MatLoad(G,viewer);
                //PetscViewerDestroy(&viewer);
                Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
                ai= aa->i;
                aj= aa->j;
                vals=aa->a;
                nz01=aa->nz;
                count=0;
                for(i=0; i<nz01; i++) if(vals[i]!=0) {
                    count++;
                  }
              }
              //MPI_Bcast(&count,1, MPI_LONG, 0, PETSC_COMM_WORLD);
              //MPI_Bcast(&nz,1, MPI_LONG, 0, PETSC_COMM_WORLD);
              int *irn=(int *) calloc (count,sizeof(int));
              int *irn1=(int *) calloc (nz01,sizeof(int));
              int *jcn=(int *) calloc (count,sizeof(int));
              ha_cgetype *values= (ha_cgetype *) calloc (count,sizeof(ha_cgetype));
              if(rank==rank_hsl) {
                for(i=0; i<VecSize-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
                    irn1[j]=i+1;
                  }
                for(j=ai[VecSize-1]; j<nz01; j++) {
                  irn1[j]=VecSize;
                }
                j=0;
                for(i=0; i<nz01; i++) if(vals[i]!=0) {
                    irn[j]=irn1[i];
                    jcn[j]=aj[i]+1;
                    values[j]=vals[i];
                    //if(i<45)printf("aj %d jcn %d\n",aj[i]+1,jcn[j]);
                    j++;
                  }
              }
              ierr = MatDestroy(&A);
              CHKERRQ(ierr);
              free(irn1);
              b1=realloc (b1,VecSize*sizeof(ha_cgetype));
              if(rank==rank_hsl) {
                VecGetArray(vecb,&vals);
                for(i=0; i<VecSize; i++) {
                  b1[i]=vals[i];
                }
              }
              ierr = VecDestroy(&vecb);
              CHKERRQ(ierr);
              int *neleperrow= (int *) calloc (VecSize,sizeof(int));
              int *ai1= (int *) calloc (VecSize,sizeof(int));
              //printf("count %d\n",count);
              if(rank==rank_hsl) {
                j=1;
                k=0,m=1;
                for(i=1; i<count; i++) {
                  if(irn[i]-irn[i-1]>0) {
                    neleperrow[k]=j;
                    ai1[k]=m;
                    j=1;
                    m=i+1;
                    k++;
                  }
                  else {
                    j++;
                  }
                }
                neleperrow[k]=j;
                ai1[k]=ai1[k-1]+neleperrow[k-1];
              }
              //ierr = MatDestroy(&A);
              //CHKERRQ(ierr);
              //free(irn1);
              //ierr = VecDestroy(&b);
              //CHKERRQ(ierr);
              //MPI_Bcast(irn, count, MPI_LONG, 0, PETSC_COMM_WORLD);
              //MPI_Bcast(jcn, count, MPI_LONG, 0, PETSC_COMM_WORLD);
              //MPI_Bcast(b1, VecSize, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
              //MPI_Bcast(values, count, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
              //MPI_Bcast(neleperrow, VecSize, MPI_LONG, 0, PETSC_COMM_WORLD);
              //MPI_Bcast(ai1, VecSize, MPI_LONG, 0, PETSC_COMM_WORLD);
              indata[1]=VecSize;//.m
              indata[0]=count;//.nz
              ptx = indata;
              x1=realloc (x1,VecSize*sizeof(ha_cgetype));
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
              //if(alltimeset>=0||allregset>=0) {
              if(mc66!=0)spec48_single_(ptx,irn,jcn,b1,values,x1,neleperrow,ai1,&fcomm);
              free(irn);
              if(mc66==0)spec48_nomc66_(ptx,jcn,b1,values,x1,neleperrow,&fcomm,counteq,countvarintra1);
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
              //MPI_Bcast(x1, VecSize, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
              //printf("rank %d\n\n",rank);
              free(jcn);
              free(values);
              free(neleperrow);
              free(ai1);
              free(b1);
              b1=NULL;
              //b1=realloc (b1,sizeof(ha_cgetype));
            }
            else {
              if(rank==rank_hsl) {
                //PetscViewerBinaryOpen(PETSC_COMM_SELF,"A1.out",FILE_MODE_READ,&viewer);
                //MatLoad(G,viewer);
                //PetscViewerDestroy(&viewer);
                Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
                ai= aa->i;
                aj= aa->j;
                vals=aa->a;
                nz01=aa->nz;
                count=0;
                for(i=0; i<nz01; i++) if(vals[i]!=0) {
                    count++;
                  }

              }
              lasize=ceil((laA/100)*count);
              int *irn=(int *) calloc (lasize,sizeof(int));
              int *irn1=(int *) calloc (nz01,sizeof(int));
              int *jcn=(int *) calloc (lasize,sizeof(int));
              ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
              //printf("count %d\n",count);
              if(rank==rank_hsl) {
                for(i=0; i<VecSize-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
                    irn1[j]=i+1;
                  }
                for(j=ai[VecSize-1]; j<nz01; j++) {
                  irn1[j]=VecSize;
                }
                j=0;
                for(i=0; i<nz01; i++) if(vals[i]!=0) {
                    irn[j]=irn1[i];
                    jcn[j]=aj[i]+1;
                    values[j]=vals[i];
                    //if(i<45)printf("aj %d jcn %d\n",aj[i]+1,jcn[j]);
                    j++;
                  }
                VecGetArray(vecb,&vals);
              }
              ierr = MatDestroy(&A);
              CHKERRQ(ierr);
              free(irn1);
              x1=realloc (x1,VecSize*sizeof(ha_cgetype));
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
              int *insize=(int *) calloc (4,sizeof(int));
              insize[0]=VecSize;
              insize[1]=VecSize;
              insize[2]=count;
              insize[3]=laA;
              //spec48_ssol_(insize,irn,jcn,values,b1,x0);
              if(rank==rank_hsl)spec48_ssol2la_(insize,irn,jcn,values,vals,x1);
              ierr = VecDestroy(&vecb);
              CHKERRQ(ierr);
              free(insize);
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
              //MPI_Bcast(x1, VecSize, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
              //printf("rank %d\n\n",rank);
              free(irn);
              free(jcn);
              free(values);
            }
          }
          //if(sol==2&&stepcount==4){
          //for(i=0;i<VecSize;i++)printf("xc %lf\n",x1[i]);
          //exit(0);
          //}
          if(rank==rank_hsl) {
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            ha_cgeshock=realloc (ha_cgeshock,(nvarele)*sizeof(ha_cgeexovar));
            freadresult=fread(ha_cgeshock, sizeof(ha_cgeexovar),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempvar");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            ha_cofvar=realloc (ha_cofvar,(ncofele+nvarele)*sizeof(ha_cgevar));
            freadresult=fread(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempclag1");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            clag1=realloc (clag1,(nvarele)*sizeof(ha_cgetype));
            freadresult=fread(clag1, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempvarchange");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            varchange=realloc (varchange,(nvarele)*sizeof(ha_cgetype));
            freadresult=fread(varchange, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);

//            strcpy(tempfilenam,temdir);
//            strcat(tempfilenam,"_tempxcf");
//            sprintf(tempchar, "%d",rank);
//            strcat(tempfilenam,tempchar);
//            strcat(tempfilenam,".bin");
//            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
//              printf("Error opening file\n");
//            }
//            xcf=realloc (xcf,(nvarele)*sizeof(ha_cgetype));
//            freadresult=fread(xcf, sizeof(ha_cgetype),nvarele, tempvar);
//            fclose(tempvar);
//            remove(tempfilenam);
//            if(sizeof(xc0)>sizeof(ha_cgetype)) {
//              strcpy(tempfilenam,temdir);
//              strcat(tempfilenam,"_tempxcO");
//              sprintf(tempchar, "%d",rank);
//              strcat(tempfilenam,tempchar);
//              strcat(tempfilenam,".bin");
//              if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
//                printf("Error opening file\n");
//              }
//              xc0=realloc (xc0,(nvarele)*sizeof(ha_cgetype));
//              freadresult=fread(xc0, sizeof(ha_cgetype),nvarele, tempvar);
//              fclose(tempvar);
//              remove(tempfilenam);
//            }
          }
          /*if(rank==rank_hsl) {
            if ((tempvar = fopen("_tempvar.bin", "rb")) == NULL) {
              printf("Error opening file\n");
            }
            ha_cofvar=realloc (ha_cofvar,(ncofele+nvarele)*sizeof(ha_cgevar));
            fread(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
            fclose(tempvar);
          }*/
          MPI_Barrier(PETSC_COMM_WORLD);
          printf("sol %d stepcount %d\n\n",sol,stepcount);
          if(nohsl) {
            VecCreate(PETSC_COMM_WORLD,&vece);
          }
          else {
            VecCreate(PETSC_COMM_SELF,&vece);
          }
          if(nohsl) {
            VecSetType(vece,VECMPI);
          }
          else {
            VecSetType(vece,VECSEQ);
          }
          printf("rank %d check %d\n\n",rank,1);
          if(nesteddbbd==1)VecSetSizes(vece,localsize,VecSize);
          else VecSetSizes(vece,PETSC_DECIDE,VecSize);
          VecSetOption(vece, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
          printf("rank %d check %d\n\n",rank,2);
          ha_cofvar1=ha_cofvar+ncofele;
          if(stepcount==0) {
            for(i=0; i<nvar; i++) {
              if(ha_var[i].change_real) {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  //tindx1=ha_var[i].begadd+j;
                  //tindx2=tindx1+ncofele;
                  if(ha_cgeshock[tindx1].ShockId) {
                    ha_cofvar1[tindx1].varval+=ha_cofvar1[tindx1].csolpupd;
                    varchange[tindx1]=ha_cofvar1[tindx1].csolpupd;
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,ha_cofvar1[tindx1].csolpupd,INSERT_VALUES);
                  }
                  else {
                    varchange[tindx1]=x1[ha_cgeshock[tindx1].ExoIndx];
                    ha_cofvar1[tindx1].varval+=x1[ha_cgeshock[tindx1].ExoIndx];
                    ha_cofvar1[tindx1].csolpupd=x1[ha_cgeshock[tindx1].ExoIndx];
                    clag1[tindx1]=0;
                  }
                }
              }
              else {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  //tindx1=ha_var[i].begadd+j;
                  //tindx2=tindx1+ncofele;
                  if(ha_cgeshock[tindx1].ShockId) {
                    varchange[tindx1]=ha_cofvar1[tindx1].csolpupd;
                    ha_cofvar1[tindx1].varval*=(1+ha_cofvar1[tindx1].csolpupd/100);
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,ha_cofvar1[tindx1].csolpupd/(1+ha_cofvar1[tindx1].csolpupd/100),INSERT_VALUES);
                  }
                  else {
                    varchange[tindx1]=x1[ha_cgeshock[tindx1].ExoIndx];
                    ha_cofvar1[tindx1].csolpupd=x1[ha_cgeshock[tindx1].ExoIndx];
                    ha_cofvar1[tindx1].varval*=(1+ha_cofvar1[tindx1].csolpupd/100);
                    clag1[tindx1]=0;
                  }
                }
              }
            }
          }
          else {
            for(i=0; i<nvar; i++) {
              if(ha_var[i].change_real) {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  //tindx1=ha_var[i].begadd+j;
                  //tindx2=tindx1+ncofele;
                  if(ha_cgeshock[tindx1].ShockId) {
                    ha_cofvar1[tindx1].varval+=ha_cofvar1[tindx1].csolpupd;
                    varchange[tindx1]+=ha_cofvar1[tindx1].csolpupd;
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,ha_cofvar1[tindx1].csolpupd,INSERT_VALUES);
                  }
                  else {
                    //temp1=ha_cofvar[ncofele+ha_var[i].begadd+j].varchange;
                    //ha_cofvar[ncofele+ha_var[i].begadd+j].varchange=clag1[ha_var[i].begadd+j]+2*x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx];//+=x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx];//
                    //ha_cofvar[ncofele+ha_var[i].begadd+j].csolpupd=x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx];//ha_cofvar[ncofele+ha_var[i].begadd+j].varchange-temp1;
                    //ha_cofvar[ncofele+ha_var[i].begadd+j].varval+=ha_cofvar[ncofele+ha_var[i].begadd+j].varchange-temp1;//ha_cofvar[ncofele+ha_var[i].begadd+j].csolpupd;
                    //clag1[ha_var[i].begadd+j]=temp1;
                    temp1=ha_cofvar1[tindx1].varval;//change;
                    varchange[tindx1]=clag1[tindx1]+2*x1[ha_cgeshock[tindx1].ExoIndx];//+=x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx];//
                    ha_cofvar1[tindx1].csolpupd=x1[ha_cgeshock[tindx1].ExoIndx];//ha_cofvar[ncofele+ha_var[i].begadd+j].varchange-temp1;
                    ha_cofvar1[tindx1].varval=clag1[tindx1]+2*x1[ha_cgeshock[tindx1].ExoIndx];//ha_cofvar[ncofele+ha_var[i].begadd+j].varchange-temp1;//ha_cofvar[ncofele+ha_var[i].begadd+j].csolpupd;
                    clag1[tindx1]=temp1;
                  }
                }
              }
              else {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  //tindx1=ha_var[i].begadd+j;
                  //tindx2=tindx1+ncofele;
                  if(ha_cgeshock[tindx1].ShockId) {
                    //temp1=ha_cgeshock[ha_var[i].begadd+j].ShockVal/nsteps;
                    temp2=ha_cgeshock[tindx1].ShockVal;//subints;
                    temp1=(100+(subindx+1)*temp2)/(100+subindx*temp2)-1;
                    temp1*=vpercents;
                    ha_cofvar1[tindx1].csolpupd=temp1/(1+varchange[tindx1]/100);
                    varchange[tindx1]+=temp1;//*(1+ha_cofvar[ncofele+ha_var[i].begadd+j].varchange/100)
                    ha_cofvar1[tindx1].varval=(1+varchange[tindx1]/100)*ha_cofvar1[tindx1].var0;
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,temp1/(1+varchange[tindx1]/100),INSERT_VALUES);
                  }
                  else {
                    //temp1=ha_cofvar[ncofele+ha_var[i].begadd+j].varchange;
                    //ha_cofvar[ncofele+ha_var[i].begadd+j].varchange=clag1[ha_var[i].begadd+j]+2*x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx]*(100+temp1)/100;//+=x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx]*(1+temp1/100);//
                    //ha_cofvar[ncofele+ha_var[i].begadd+j].csolpupd=x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx];//(ha_cofvar[ncofele+ha_var[i].begadd+j].varchange-temp1)/(1+temp1/100);
                    //ha_cofvar[ncofele+ha_var[i].begadd+j].varval=(100+ha_cofvar[ncofele+ha_var[i].begadd+j].varchange)/100*ha_cofvar[ncofele+ha_var[i].begadd+j].var0;
                    //clag1[ha_var[i].begadd+j]=temp1;
                    temp1=varchange[tindx1];
                    varchange[tindx1]=clag1[tindx1]+2*x1[ha_cgeshock[tindx1].ExoIndx]*(100+temp1)/100;//+=x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx]*(1+temp1/100);//
                    ha_cofvar1[tindx1].csolpupd=x1[ha_cgeshock[tindx1].ExoIndx];//(ha_cofvar[ncofele+ha_var[i].begadd+j].varchange-temp1)/(1+temp1/100);
                    ha_cofvar1[tindx1].varval=(100+varchange[tindx1])/100*ha_cofvar1[tindx1].var0;
                    clag1[tindx1]=temp1;
                  }
                }
              }
            }
          }
          //printf("X1\n");
          printf("rank %d check %d\n\n",rank,3);
          free(x1);
          x1=NULL;
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = VecAssemblyBegin(vece);
          CHKERRQ(ierr);
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = VecAssemblyEnd(vece);
          CHKERRQ(ierr);
          //x1=realloc (x1,sizeof(ha_cgetype));
          printf("rank %d check %d\n\n",rank,4);
          if(rank==rank_hsl) {
            if(stepcount==0) {
              hnew_update(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele);
            }
            else {
              hnew_mupdate(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele);
            }
            //printf("X2\n");
            strcpy(commsyntax,"formula");
            IsIni=false;
            hnew_calcff(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,IsIni);
            //printf("X3\n");
          }
          printf("rank %d check %d\n\n",rank,5);
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time1);
          CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Update time %f\n",time1-time0);
          CHKERRQ(ierr);
          ierr = PetscGetCPUTime(&time0);
          CHKERRQ(ierr);
          printf("OKKL!\n");
        }

        strcpy(commsyntax,"equation");
        //ierr = MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize,dnz,dnnz,onz,onnz,&A);
        //CHKERRQ(ierr);
        //ierr = MatSetFromOptions(A);
        //CHKERRQ(ierr);
        //ierr = MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize,dnzB,dnnzB,onzB,onnzB,&B);
        //CHKERRQ(ierr);
        //ierr = MatSetFromOptions(B);
        //CHKERRQ(ierr);
        if(rank==rank_hsl) {
//          if(sizeof(xc0)>sizeof(ha_cgetype)) {
//            strcpy(tempfilenam,temdir);
//            strcat(tempfilenam,"_tempxcO");
//            sprintf(tempchar, "%d",rank);
//            strcat(tempfilenam,tempchar);
//            strcat(tempfilenam,".bin");
//            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
//              printf("Error opening file\n");
//              //return 1;
//            }
//            fwrite(xc0, sizeof(ha_cgetype),nvarele, tempvar);
//            fclose(tempvar);
//            xc0=realloc (xc0,2*sizeof(ha_cgetype));
//          }

//          strcpy(tempfilenam,temdir);
//          strcat(tempfilenam,"_tempxcf");
//          sprintf(tempchar, "%d",rank);
//          strcat(tempfilenam,tempchar);
//          strcat(tempfilenam,".bin");
//          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
//            printf("Error opening file\n");
//            //return 1;
//          }
//          fwrite(xcf, sizeof(ha_cgetype),nvarele, tempvar);
//          fclose(tempvar);
//          xcf=realloc (xcf,1*sizeof(ha_cgetype));

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempclag1");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
            printf("Error opening file\n");
            //return 1;
          }
          fwrite(clag1, sizeof(ha_cgetype),nvarele, tempvar);
          fclose(tempvar);
          free(clag1);
          clag1=NULL;
          //clag1=realloc (clag1,1*sizeof(ha_cgetype));

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempvarchange");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
            printf("Error opening file\n");
            //return 1;
          }
          fwrite(varchange, sizeof(ha_cgetype),nvarele, tempvar);
          fclose(tempvar);
          free(varchange);
          varchange=NULL;
          //varchange=realloc (varchange,1*sizeof(ha_cgetype));

        }

        MPI_Barrier(PETSC_COMM_WORLD);
        if(nohsl) {
          MatCreate(PETSC_COMM_WORLD,&A);
        }
        else {
          MatCreate(PETSC_COMM_SELF,&A);
        }
        if(nesteddbbd==1)MatSetSizes(A,localsize,localsize,VecSize,VecSize);
        else MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
        if(nohsl) {
          MatSetType(A,MATMPIAIJ);
        }
        else {
          MatSetType(A,MATSEQAIJ);
        }
        if(nohsl) {
          MatMPIAIJSetPreallocation(A,dnz,dnnz,onz,onnz);
        }
        else {
          MatSeqAIJSetPreallocation(A,dnz,dnnz);
        }

        if(nohsl) {
          MatCreate(PETSC_COMM_WORLD,&B);
        }
        else {
          MatCreate(PETSC_COMM_SELF,&B);
        }
        if(nesteddbbd==1)MatSetSizes(B,localsize,localsize,VecSize,VecSize);
        else MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
        if(nohsl) {
          MatSetType(B,MATMPIAIJ);
        }
        else {
          MatSetType(B,MATSEQAIJ);
        }
        if(nohsl) {
          MatMPIAIJSetPreallocation(B,dnzB,dnnzB,onzB,onnzB);
        }
        else {
          MatSeqAIJSetPreallocation(B,dnzB,dnnzB);
        }

        if(rank==rank_hsl) {
          HaNewMatVal(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,ndblock,alltimeset,allregset,ha_eqadd,counteq,nintraeq,A,B);
        }

        if(rank==rank_hsl) {
          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempvar");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
            printf("Error opening file\n");
            //return 1;
          }
          fwrite(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
          fclose(tempvar);
          free(ha_cofvar);
          ha_cofvar=NULL;
          //ha_cofvar=realloc (ha_cofvar,1*sizeof(ha_cgevar));

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempshock");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
            printf("Error opening file\n");
            //return 1;
          }
          fwrite(ha_cgeshock, sizeof(ha_cgeexovar),nvarele, tempvar);
          fclose(tempvar);
          free(ha_cgeshock);
          ha_cgeshock=NULL;
          //ha_cgeshock=realloc (ha_cgeshock,1*sizeof(ha_cgeexovar));
        }

        MPI_Barrier(PETSC_COMM_WORLD);
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        ierr = VecDuplicate(vece,&vecb);
        CHKERRQ(ierr);
        if(rank==rank_hsl) {
          ierr = MatMult(B,vece,vecb);
          CHKERRQ(ierr);
        }
        ierr = VecDestroy(&vece);
        CHKERRQ(ierr);
        ierr = VecAssemblyBegin(vecb);
        CHKERRQ(ierr);
        ierr = VecAssemblyEnd(vecb);
        CHKERRQ(ierr);
        /*PetscViewerBinaryOpen(PETSC_COMM_WORLD, "b.out",FILE_MODE_WRITE,&viewer);
        ierr = VecView(b,viewer);
        CHKERRQ(ierr);
        PetscViewerDestroy(&viewer);
        VecCreateSeq(PETSC_COMM_SELF,VecSize,&b2);
        if(rank==0) {
          PetscViewerBinaryOpen(PETSC_COMM_SELF,"b.out",FILE_MODE_READ,&viewer);
          VecLoad(b2,viewer);
          PetscViewerDestroy(&viewer);
        }
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, "A1.out",FILE_MODE_WRITE, &viewer);
        ierr = MatView(A,viewer);
        CHKERRQ(ierr);
        PetscViewerDestroy(&viewer);
        ierr = MatDestroy(&A);
        CHKERRQ(ierr);*/
        ierr = MatDestroy(&B);
        CHKERRQ(ierr);
        //MatCreate(PETSC_COMM_SELF,&G);
        //MatSetType(G,MATSEQAIJ);
        /*if(rank==rank_hsl) {
          if ( (tempvar = fopen("_tempvar.bin", "wb")) == NULL ) {
            printf("Error opening file\n");
            //return 1;
          }
          fwrite(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
          fclose(tempvar);
          ha_cofvar=realloc (ha_cofvar,sizeof(ha_cgevar));
        }*/
        if(matsol>=2) {
          //strcpy(tempchar,"lrzip -D -l -L 1 ");
          //strcat(tempchar,tempfilenam);
          //system(tempchar);
          //remove(tempfilenam);
          //ha_cofvar=realloc (ha_cofvar,sizeof(ha_cgevar));
          int *ha_rows= (int *) calloc (VecSize,sizeof(int));
          int *ha_cols= (int *) calloc (VecSize,sizeof(int));
          int *ha_ndblocks= (int *) calloc (ndblock,sizeof(int));
          //ha_eqadd1=realloc(ha_eqadd1,VecSize*sizeof(uvadd));
          //memcpy (ha_eqadd1,ha_eqadd,VecSize*sizeof(uvadd));
          //uvadd *orig_exoindx= (uvadd *) calloc (VecSize,sizeof(uvadd));
          //timestr=clock();//time(&timestr);
          gettimeofday(&gettime_now, NULL);

          if(matsol==2) {
            HaDBBDMatOder(A,VecSize,mpisize,rank,Istart,Iend,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,cntl6);
            x1=realloc (x1,VecSize*sizeof(ha_cgetype));
            //VecView(b,0);
            //VecSet(b,1);
            HaDBBDParSol(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laD,cntl3);//,iter
          }

          if(matsol==3) {
            presol=1;
//            uvadd *counteqs= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
//            uvadd *counteqnoadds= (uvadd *) calloc (ndblock,sizeof(uvadd));
//            uvadd *countvarintra1s= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
            memcpy(counteq,counteqs,(ndblock+1)*sizeof(uvadd));
            memcpy(counteqnoadd,counteqnoadds,(ndblock)*sizeof(uvadd));
            memcpy(countvarintra1,countvarintra1s,(ndblock+1)*sizeof(uvadd));
            HaNDBBDMatOderPre(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
            HaNDBBDParPre(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,nreg,ntime,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
//        printf("Would you like to record ranks?\n");
//        ch = getchar();
            presol=0;
//            memcpy(counteqs,counteq,(ndblock+1)*sizeof(uvadd));
//            memcpy(counteqnoadds,counteqnoadd,(ndblock)*sizeof(uvadd));
//            memcpy(countvarintra1s,countvarintra1,(ndblock+1)*sizeof(uvadd));
            HaNDBBDMatOder(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
            x1=realloc (x1,VecSize*sizeof(ha_cgetype));
            HaNDBBDParSol(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,nreg,ntime,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
//            free(counteqs);
//            free(counteqnoadds);
//            free(countvarintra1s);
          }

          //timeend=clock();//time(&timeend);
          gettimeofday(&endtime, NULL);
          //memcpy (ha_eqadd,ha_eqadd1,VecSize*sizeof(uvadd));
          //ha_eqadd1=realloc(ha_eqadd1,1*sizeof(uvadd));
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time1);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"One step solution %f\n",time1-time0);
          //if(rank==0)printf("One step calculation time %f\n",difftime(timeend,timestr)/CLOCKS_PER_SEC);
          if(rank==0)printf("One step calculation time %f\n",(endtime.tv_sec - gettime_now.tv_sec)+((double)(endtime.tv_usec - gettime_now.tv_usec))/ 1000000);
          free(ha_rows);
          free(ha_cols);
          free(ha_ndblocks);
          printf("rank %d\n",rank);
          MPI_Barrier(PETSC_COMM_WORLD);
        }
        else {
          if(matsol==1) {
            if(rank==rank_hsl) {
              //PetscViewerBinaryOpen(PETSC_COMM_SELF,"A1.out",FILE_MODE_READ,&viewer);
              //MatLoad(G,viewer);
              //PetscViewerDestroy(&viewer);
              Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
              ai= aa->i;
              aj= aa->j;
              vals=aa->a;
              nz01=aa->nz;
              count=0;
              for(i=0; i<nz01; i++) if(vals[i]!=0) {
                  count++;
                }
            }
            //MPI_Bcast(&count,1, MPI_LONG, 0, PETSC_COMM_WORLD);
            //MPI_Bcast(&nz,1, MPI_LONG, 0, PETSC_COMM_WORLD);
            int *irn=(int *) calloc (count,sizeof(int));
            int *irn1=(int *) calloc (nz01,sizeof(int));
            int *jcn=(int *) calloc (count,sizeof(int));
            ha_cgetype *values= (ha_cgetype *) calloc (count,sizeof(ha_cgetype));
            if(rank==rank_hsl) {
              for(i=0; i<VecSize-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
                  irn1[j]=i+1;
                }
              for(j=ai[VecSize-1]; j<nz01; j++) {
                irn1[j]=VecSize;
              }
              j=0;
              for(i=0; i<nz01; i++) if(vals[i]!=0) {
                  irn[j]=irn1[i];
                  jcn[j]=aj[i]+1;
                  values[j]=vals[i];
                  //if(i<45)printf("aj %d jcn %d\n",aj[i]+1,jcn[j]);
                  j++;
                }
            }
            ierr = MatDestroy(&A);
            CHKERRQ(ierr);
            free(irn1);
            b1=realloc (b1,VecSize*sizeof(ha_cgetype));
            if(rank==rank_hsl) {
              VecGetArray(vecb,&vals);
              for(i=0; i<VecSize; i++) {
                b1[i]=vals[i];
              }
            }
            ierr = VecDestroy(&vecb);
            CHKERRQ(ierr);
            int *neleperrow= (int *) calloc (VecSize,sizeof(int));
            int *ai1= (int *) calloc (VecSize,sizeof(int));
            if(rank==rank_hsl) {
              j=1;
              k=0,m=1;
              for(i=1; i<count; i++) {
                if(irn[i]-irn[i-1]>0) {
                  neleperrow[k]=j;
                  ai1[k]=m;
                  j=1;
                  m=i+1;
                  k++;
                }
                else {
                  j++;
                }
              }
              neleperrow[k]=j;
              ai1[k]=ai1[k-1]+neleperrow[k-1];
              //for(i=0;i<VecSize;i++)printf("nele %d ai1 %d row %d row1 %d\n",neleperrow[i],ai1[i],i,irn[ai1[i]]);
              //for(i=count-100; i<count; i++) printf("x %f irn %d jrn %d count %d nz %d\n",values[i],irn[i],jcn[i],count,nz);
            }
            //ierr = VecDestroy(&b);
            //CHKERRQ(ierr);
            //MPI_Bcast(irn, count, MPI_LONG, 0, PETSC_COMM_WORLD);
            //MPI_Bcast(jcn, count, MPI_LONG, 0, PETSC_COMM_WORLD);
            //MPI_Bcast(b1, VecSize, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
            //MPI_Bcast(values, count, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
            //MPI_Bcast(neleperrow, VecSize, MPI_LONG, 0, PETSC_COMM_WORLD);
            //MPI_Bcast(ai1, VecSize, MPI_LONG, 0, PETSC_COMM_WORLD);
            indata[1]=VecSize;//.m
            indata[0]=count;//.nz
            ptx = indata;//&
            x1=realloc (x1,VecSize*sizeof(ha_cgetype));
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
            //if(alltimeset>=0||allregset>=0) {
            if(mc66!=0)spec48_single_(ptx,irn,jcn,b1,values,x1,neleperrow,ai1,&fcomm);
            free(irn);
            if(mc66==0)spec48_nomc66_(ptx,jcn,b1,values,x1,neleperrow,&fcomm,counteq,countvarintra1);
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
            //MPI_Bcast(x1, VecSize, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
            free(jcn);
            free(values);
            free(neleperrow);
            free(ai1);
            free(b1);
            b1=NULL;
            //b1=realloc (b1,sizeof(ha_cgetype));
          }
          else {
            if(rank==rank_hsl) {
              //PetscViewerBinaryOpen(PETSC_COMM_SELF,"A1.out",FILE_MODE_READ,&viewer);
              //MatLoad(G,viewer);
              //PetscViewerDestroy(&viewer);
              Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
              ai= aa->i;
              aj= aa->j;
              vals=aa->a;
              nz01=aa->nz;
              count=0;
              for(i=0; i<nz01; i++) if(vals[i]!=0) {
                  count++;
                }
            }
            //MPI_Bcast(&count,1, MPI_LONG, 0, PETSC_COMM_WORLD);
            //MPI_Bcast(&nz,1, MPI_LONG, 0, PETSC_COMM_WORLD);
            lasize=ceil((laA/100)*count);
            int *irn=(int *) calloc (lasize,sizeof(int));
            int *irn1=(int *) calloc (nz01,sizeof(int));
            int *jcn=(int *) calloc (lasize,sizeof(int));
            ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
            if(rank==rank_hsl) {
              for(i=0; i<VecSize-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
                  irn1[j]=i+1;
                }
              for(j=ai[VecSize-1]; j<nz01; j++) {
                irn1[j]=VecSize;
              }
              j=0;
              for(i=0; i<nz01; i++) if(vals[i]!=0) {
                  irn[j]=irn1[i];
                  jcn[j]=aj[i]+1;
                  values[j]=vals[i];
                  //if(i<45)printf("aj %d jcn %d\n",aj[i]+1,jcn[j]);
                  j++;
                }
              VecGetArray(vecb,&vals);
            }
            ierr = MatDestroy(&A);
            CHKERRQ(ierr);
            free(irn1);
            x1=realloc (x1,VecSize*sizeof(ha_cgetype));
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
            int *insize=(int *) calloc (4,sizeof(int));
            insize[0]=VecSize;
            insize[1]=VecSize;
            insize[2]=count;
            insize[3]=laA;
            //spec48_ssol_(insize,irn,jcn,values,b1,x0);
            if(rank==rank_hsl)spec48_ssol2la_(insize,irn,jcn,values,vals,x1);
            free(insize);
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
            //MPI_Bcast(x1, VecSize, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
            ierr = VecDestroy(&vecb);
            CHKERRQ(ierr);
            free(irn);
            free(jcn);
            free(values);
          }
        }
        if(rank==rank_hsl) {
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error opening file\n");
          }
          ha_cgeshock=realloc (ha_cgeshock,(nvarele)*sizeof(ha_cgeexovar));
          freadresult=fread(ha_cgeshock, sizeof(ha_cgeexovar),nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempvar");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error opening file\n");
          }
          ha_cofvar=realloc (ha_cofvar,(ncofele+nvarele)*sizeof(ha_cgevar));
          freadresult=fread(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempclag1");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error opening file\n");
          }
          clag1=realloc (clag1,(nvarele)*sizeof(ha_cgetype));
          freadresult=fread(clag1, sizeof(ha_cgetype),nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempvarchange");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error opening file\n");
          }
          varchange=realloc (varchange,(nvarele)*sizeof(ha_cgetype));
          freadresult=fread(varchange, sizeof(ha_cgetype),nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);

//          strcpy(tempfilenam,temdir);
//          strcat(tempfilenam,"_tempxcf");
//          sprintf(tempchar, "%d",rank);
//          strcat(tempfilenam,tempchar);
//          strcat(tempfilenam,".bin");
//          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
//            printf("Error opening file\n");
//          }
//          xcf=realloc (xcf,(nvarele)*sizeof(ha_cgetype));
//          freadresult=fread(xcf, sizeof(ha_cgetype),nvarele, tempvar);
//          fclose(tempvar);
//          remove(tempfilenam);
//          if(sizeof(xc0)>sizeof(ha_cgetype)) {
//            strcpy(tempfilenam,temdir);
//            strcat(tempfilenam,"_tempxcO");
//            sprintf(tempchar, "%d",rank);
//            strcat(tempfilenam,tempchar);
//            strcat(tempfilenam,".bin");
//            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
//              printf("Error opening file\n");
//            }
//            xc0=realloc (xc0,(nvarele)*sizeof(ha_cgetype));
//            freadresult=fread(xc0, sizeof(ha_cgetype),nvarele, tempvar);
//            fclose(tempvar);
//            remove(tempfilenam);
//          }
        }
        /*if(rank==rank_hsl) {
          if ((tempvar = fopen("_tempvar.bin", "rb")) == NULL) {
            printf("Error opening file\n");
          }
          ha_cofvar=realloc (ha_cofvar,(ncofele+nvarele)*sizeof(ha_cgevar));
          fread(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
          fclose(tempvar);
        }*/
        ha_cofvar1=ha_cofvar+ncofele;
        for(i=0; i<nvar; i++) {
          if(ha_var[i].change_real) {
            for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
              //tindx1=ha_var[i].begadd+j;
              //tindx2=tindx1+ncofele;
              if(ha_cgeshock[tindx1].ShockId) {
                ha_cofvar1[tindx1].varval=0;
                //ha_cofvar1[tindx1].csolpupd=0;
                //varchange[tindx1]+=ha_cofvar[tindx2].csolpupd;
              }
              else {
                varchange[tindx1]=0.5*(varchange[tindx1]+clag1[tindx1]+x1[ha_cgeshock[tindx1].ExoIndx]);
                ha_cofvar1[tindx1].varval=0;//ha_cofvar[tindx2].var0+varchange[tindx1];//no distortion between steps
                //ha_cofvar1[tindx1].csolpupd=0;//x1[ha_cgeshock[tindx1].ExoIndx];
                clag1[tindx1]=0;
//                temp1=ha_cofvar[tindx2].varval;//change;
//                varchange[tindx1]=clag1[tindx1]+2*x1[ha_cgeshock[tindx1].ExoIndx];
//                ha_cofvar[tindx2].varval=clag1[tindx1]+2*x1[ha_cgeshock[tindx1].ExoIndx];
//                ha_cofvar[tindx2].csolpupd=x1[ha_cgeshock[tindx1].ExoIndx];
//                clag1[tindx1]=temp1;
              }
            }
          }
          else {
            for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
              //tindx1=ha_var[i].begadd+j;
              //tindx2=tindx1+ncofele;
              if(ha_cgeshock[tindx1].ShockId) {
                //temp2=ha_cgeshock[tindx1].ShockVal;//subints;
                //temp1=(100+(subindx+1)*temp2)/(100+subindx*temp2)-1;
                //temp1*=vpercents;
                //ha_cofvar1[tindx1].csolpupd=0;
                //varchange[tindx1]+=temp1;
                ha_cofvar1[tindx1].varval=0;
              }
              else {
//                temp1=varchange[tindx1];
//                varchange[tindx1]=clag1[tindx1]+2*x1[ha_cgeshock[tindx1].ExoIndx]*(100+temp1)/100;
//                ha_cofvar[tindx2].csolpupd=x1[ha_cgeshock[tindx1].ExoIndx];
//                ha_cofvar[tindx2].varval=(100+varchange[tindx1])/100*ha_cofvar[tindx2].var0;
//                clag1[tindx1]=temp1;
                varchange[tindx1]=0.5*(varchange[tindx1]+clag1[tindx1]+x1[ha_cgeshock[tindx1].ExoIndx]*(1+varchange[tindx1]/100));
                ha_cofvar1[tindx1].varval=0;//ha_cofvar[tindx2].varval*varchange[tindx1]/100;
                //ha_cofvar1[tindx1].csolpupd=0;//x1[ha_cgeshock[tindx1].ExoIndx];
                clag1[tindx1]=0;
              }
            }
          }
        }
//        x1=realloc (x1,sizeof(ha_cgetype));
        if(rank==rank_hsl) {
//           hnew_graggupd(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele);
//           strcpy(commsyntax,"formula");
//           IsIni=false;
//           hnew_calcff(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,IsIni);
// 
//           ierr = PetscGetCPUTime(&time1);
//           CHKERRQ(ierr);
//           ierr = PetscPrintf(PETSC_COMM_WORLD,"Update time %f\n",time1-time0);
//           CHKERRQ(ierr);
          /*if(sol==0)for(i=0; i<ncofele+nvarele; i++) {
              xs[i]=ha_cofvar[i].varval/45;
            }
          if(sol==1)for(i=0; i<ncofele+nvarele; i++) {
              xs[i]-=20*ha_cofvar[i].varval/45;
            }
          if(sol==2)for(i=0; i<ncofele+nvarele; i++) {
              xs[i]+=64*ha_cofvar[i].varval/45;
            }*/
          //printf("HRE!!!!!!!!!!!!!!!!!!!!!!!!!! sol %d step %d\n",sol,stepcount);
          if(subindx!=0||sol!=0) {
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxcf");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            xcf=realloc (xcf,(nvarele)*sizeof(ha_cgetype));
            freadresult=fread(xcf, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc12");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            xc12=realloc (xc12,(nvarele)*sizeof(ha_cgetype));
            freadresult=fread(xc12, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc24");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            xc24=realloc (xc24,(nvarele)*sizeof(ha_cgetype));
            freadresult=fread(xc24, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            
            xc0=realloc (xc0,(nvarele)*sizeof(ha_cgetype));
            if(subindx>0&&sol>0){
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxcO");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            freadresult=fread(xc0, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            if(subindx==subints-1&&sol==maxsol-1)remove(tempfilenam);
            }

//             xc012=realloc (xc012,(nvarele)*sizeof(ha_cgetype));
//             if(subindx>0&&sol>0){
//             strcpy(tempfilenam,temdir);
//             strcat(tempfilenam,"_tempxcO12");
//             sprintf(tempchar, "%d",rank);
//             strcat(tempfilenam,tempchar);
//             strcat(tempfilenam,".bin");
//             if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
//               printf("Error opening file\n");
//             }
//             freadresult=fread(xc012, sizeof(ha_cgetype),nvarele, tempvar);
//             fclose(tempvar);
//             if(subindx==subints-1&&sol==maxsol-1)remove(tempfilenam);
//             xc024=realloc (xc024,(nvarele)*sizeof(ha_cgetype));
//             }
// 
//             if(subindx>0&&sol>1){
//             strcpy(tempfilenam,temdir);
//             strcat(tempfilenam,"_tempxcO24");
//             sprintf(tempchar, "%d",rank);
//             strcat(tempfilenam,tempchar);
//             strcat(tempfilenam,".bin");
//             if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
//               printf("Error opening file\n");
//             }
//             freadresult=fread(xc024, sizeof(ha_cgetype),nvarele, tempvar);
//             fclose(tempvar);
//             if(subindx==subints-1&&sol==maxsol-1)remove(tempfilenam);
//             }
          }
          if(subindx==0&&sol==0) {
            xc0=realloc (xc0,nvarele*sizeof(ha_cgetype));
            xcf=realloc (xcf,nvarele*sizeof(ha_cgetype));
            //xc12=realloc (xc12,nvarele*sizeof(ha_cgetype));
            for(i=0; i<nvarele; i++)xcf[i]=0;
          }
          //if(subindx>0&&sol==0)xc012=realloc (xc012,nvarele*sizeof(ha_cgetype));
          //if(subindx>0&&sol==1)xc024=realloc (xc024,nvarele*sizeof(ha_cgetype));
          if(sol==0)xc12=realloc (xc12,nvarele*sizeof(ha_cgetype));
          if(sol==0)xc24=realloc (xc24,nvarele*sizeof(ha_cgetype));
          if(subindx>0) {
            //if(subindx==1&&sol==0)xc0=realloc (xc0,nvarele*sizeof(ha_cgetype));
            //printf("HRE sol %d step %d\n",sol,stepcount);
            if(sol==0)for(i=0; i<nvarele; i++) xc0[i]=1+xcf[i]/100;//if(i==1287)printf("sol!!!!!!!!!!!!!!!!!! %d step %d xc %lf xc0 %lf k %d\n",sol,stepcount,1.0+xc[k]/100,xc0[i],i);}
//             if(sol==0)for(i=0; i<nvarele; i++) xc012[i]=1+xc12[i]/100;
//             if(sol==1)for(i=0; i<nvarele; i++) xc024[i]=1+xc24[i]/100;
            if(sol==0) {
              kval1=1.0/(kindx1*kindx1-1.0);
              kval2=1.0/(1-kindx1*kindx1)/(1.0-kindx2*kindx2);
              for(i=0; i<nvar; i++) {
                if(ha_var[i].change_real) {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    //if(ha_cgeshock[ha_var[i].begadd+j].ShockId==0) {
                    //k=ha_var[i].begadd+j;
                    //printf("k %ld xcf %lf change %lf var %lf\n",k,xcf[k],varchange[k],ha_cofvar[k+ncofele].varval);
                    //xc12[k]=xcf[k]-varchange[k]/3;
                    xc12[k]=xcf[k]-varchange[k]*kval1;
                    xc24[k]=xcf[k];
                    //xcf[k]+=varchange[k]/45;
                    xcf[k]+=varchange[k]*kval2;
                    //}
                  }
                }
                else {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    //if(ha_cgeshock[ha_var[i].begadd+j].ShockId==0) {
                    //k=ha_var[i].begadd+j;
                    //if(k==1287)printf("sol!!!!!!!!!!!!!!!!!! %d step %d xc %lf xc0 %lf k %d\n",sol,stepcount,xc[k],xc0[k],k);
                    xc24[k]=xcf[k];
                    //xc12[k]=xcf[k]-varchange[k]*xc0[k]/3;
                    xc12[k]=xcf[k]-varchange[k]*xc0[k]*kval1;
                    //xcf[k]+=varchange[k]/45*xc0[k];
                    xcf[k]+=varchange[k]*kval2*xc0[k];//(100+xc0[k])*(100+varchange[k]/45)/100-100;//varchange[k]/45;
                    //if(ha_cgeshock[ha_var[i].begadd+j].ShockId==1&&xc[k]>0)printf("xc %lf xc0 %lf k %d\n",xc[k],xc0[k],k);
                    //}
                  }
                }
              }
            }
            if(sol==1) {
              kval1=kindx1*kindx1/(kindx1*kindx1-1.0);
              kval2=kindx1*kindx1/(kindx2*kindx2-kindx1*kindx1);
              kval3=kindx1*kindx1*kindx1*kindx1/(kindx1*kindx1-kindx2*kindx2)/(1.0-kindx1*kindx1);
              for(i=0; i<nvar; i++) {
                if(ha_var[i].change_real) {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    //k=ha_var[i].begadd+j;
                    //xc24[k]-=varchange[k]/3;
                    //xc12[k]+=4*varchange[k]/3;
                    //xcf[k]-=20*varchange[k]/45;
                    xc24[k]-=varchange[k]*kval2;
                    xc12[k]+=varchange[k]*kval1;
                    xcf[k]-=varchange[k]*kval3;
                  }
                }
                else {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    //k=ha_var[i].begadd+j;
                    //xc24[k]-=varchange[k]*xc0[k]/3;
                    //xc12[k]+=4*varchange[k]*xc0[k]/3;
                    //xcf[k]-=20*varchange[k]/45*xc0[k];
                    xc24[k]-=varchange[k]*xc0[k]*kval2;
                    xc12[k]+=varchange[k]*xc0[k]*kval1;
                    xcf[k]-=varchange[k]*xc0[k]*kval3;//(100+xc0[k])*(100-20*varchange[k]/45)/100-100;
                    //xc[k]-=20*varchange[k]/45;
                  }
                }
              }
            }
            if(sol==2) {
              kval2=kindx2*kindx2/(kindx2*kindx2-kindx1*kindx1);
              kval3=kindx2*kindx2*kindx2*kindx2/(kindx1*kindx1-kindx2*kindx2)/(1.0-kindx2*kindx2);
              for(i=0; i<nvar; i++) {
                if(ha_var[i].change_real) {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    //k=ha_var[i].begadd+j;
                    //xc24[k]+=4*varchange[k]/3;
                    //xcf[k]+=64*varchange[k]/45;
                    xc24[k]+=varchange[k]*kval2;
                    xcf[k]+=varchange[k]*kval3;
                    //ha_cofvar1[k].varval=xcf[k];
                  }
                }
                else {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    //k=ha_var[i].begadd+j;
                    //xc24[k]+=4*varchange[k]*xc0[k]/3;
                    //xcf[k]+=64*varchange[k]/45*xc0[k];
                    xc24[k]+=varchange[k]*xc0[k]*kval2;
                    xcf[k]+=varchange[k]*xc0[k]*kval3;//(100+xc0[k])*(100+64*varchange[k]/45)/100-100;
                    //ha_cofvar1[k].varval=xcf[k]*xc0[k]/100;
                    //ha_cofvar1[k].csolpupd=xcf[k];
                    //xc[k]+=64*varchange[k]/45;
                  }
                }
              }
            }
            //if(subindx==subints-1&&sol==2&&stepcount==nsteps-1)xc0=realloc (xc0,sizeof(ha_cgetype));
          }
          else {
            if(sol==0){
              kval1=1.0/(kindx1*kindx1-1);
              kval2=1.0/(1-kindx1*kindx1)/(1-kindx2*kindx2);
              for(i=0; i<nvarele; i++) {
                //printf("i %ld xc %lf\n",i,xc[i]);
                //printf("xcf %lf change %lf var %lf\n",xcf[i],varchange[i],ha_cofvar[i+ncofele].varval);
                //xc12[i]=-varchange[i]/3;
                //xcf[i]+=varchange[i]/45;
                xc12[i]=-varchange[i]*kval1;
                xcf[i]+=varchange[i]*kval2;
              }
            }
            if(sol==1) {
              kval1=kindx1*kindx1/(kindx1*kindx1-1.0);
              kval2=kindx1*kindx1/(kindx2*kindx2-kindx1*kindx1);
              kval3=kindx1*kindx1*kindx1*kindx1/(kindx1*kindx1-kindx2*kindx2)/(1.0-kindx1*kindx1);
//              if(subindx<10)strcpy(solchar,"_tempsol0");
//              else strcpy(solchar,"_tempsol");
//              sprintf(tempchar, "%d", subindx);
//              strcat(solchar,tempchar);
//              sprintf(tempchar, "%d", sol-1);
//              strcat(solchar,tempchar);
//              strcat(solchar,".bin");
//              printf("solchar %s\n",solchar);
//              if ( (solution = fopen(solchar, "wb")) == NULL ) {
//                printf("Error opening file\n");
//                return 1;
//              }
//              freadresult=fread(x1, sizeof(ha_cgetype),nvarele, solution);
//              fclose(solution);
              for(i=0; i<nvarele; i++) {
                //printf("i %ld xc %lf\n",i,xc[i]);
                //xc24[i]=-varchange[i]/3;
                //xc12[i]+=4*varchange[i]/3;
                //xcf[i]-=20*varchange[i]/45;
                xc24[i]=-varchange[i]*kval2;
                xc12[i]+=varchange[i]*kval1;
                xcf[i]-=varchange[i]*kval3;
//                x1[i]=varchange[i]-x1[i];
              }
//              strcpy(solchar,"_tempxacfile");
//              if ( (solution = fopen(solchar, "wb")) == NULL ) {
//                printf("Error opening file\n");
//                return 1;
//              }
//              fwrite(x1, sizeof(ha_cgetype),nvarele, solution);
//              fclose(solution);
            }
            if(sol==2) {
              kval2=kindx2*kindx2/(kindx2*kindx2-kindx1*kindx1);
              kval3=kindx2*kindx2*kindx2*kindx2/(kindx1*kindx1-kindx2*kindx2)/(1.0-kindx2*kindx2);
//              if(subindx<10)strcpy(solchar,"_tempsol0");
//              else strcpy(solchar,"_tempsol");
//              sprintf(tempchar, "%d", subindx);
//              strcat(solchar,tempchar);
//              sprintf(tempchar, "%d", sol-2);
//              strcat(solchar,tempchar);
//              strcat(solchar,".bin");
//              printf("solchar %s\n",solchar);
//              if ( (solution = fopen(solchar, "wb")) == NULL ) {
//                printf("Error opening file\n");
//                return 1;
//              }
//              freadresult=fread(x1, sizeof(ha_cgetype),nvarele, solution);
//              fclose(solution);
              for(i=0; i<nvarele; i++) {
                //printf("i %ld xc %lf var %lf\n",i,varchange[i],xc[i]);
                //xc24[i]+=4*varchange[i]/3;
                //xcf[i]+=64*varchange[i]/45;
                xc24[i]+=varchange[i]*kval2;
                xcf[i]+=varchange[i]*kval3;
                //ha_cofvar1[i].varval=xcf[i];
                //ha_cofvar1[i].csolpupd=xcf[i];
                //printf("i %ld xc %lf var %lf\n",i,varchange[i],xc[i]);
//                x1[i]=varchange[i]-x1[i];
              }
            }
          }

          if(sol==maxsol-1){
          if(subindx==0){
          for(i=0; i<nvar; i++) {
            //printf("i %ld var %s indx %ld\n",i,ha_var[i].cofname,ncofele+ha_var[i].begadd);
            for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
               //ha_cofvar1[tindx1].varval=100;//xcf[ha_cgeshock[tindx1].ExoIndx];
               ha_cofvar1[tindx1].csolpupd=xcf[tindx1];
               //if(rank==0&&i==0)printf("i %ld var %lf xsf %lf\n",tindx1,ha_cofvar1[tindx1].csolpupd,xcf[tindx1]);
            }
          }
          }else{
          for(i=0; i<nvar; i++) {
            for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
               ha_cofvar1[tindx1].csolpupd=(100+xcf[tindx1])/xc0[tindx1]-100;
            }
          }
          }
//           for(i=ncofele; i<ncofele+nvarele; i++){
//             ha_cofvar[i].csolpupd=xcf[i];
//             if(i<ncofele+6)printf("i %d var %lf\n",i,ha_cofvar[i].varval);
//           }
          for(i=0; i<ncofele; i++) ha_cofvar[i].varval=ha_cofvar[i].var0;
          hnew_gupd(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele);
          strcpy(commsyntax,"formula");
          IsIni=false;
          hnew_calcff(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,IsIni);
          for(i=0; i<nvar; i++) {
            //printf("i %ld var %s indx %ld\n",i,ha_var[i].cofname,ncofele+ha_var[i].begadd);
            for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
               //ha_cofvar1[tindx1].varval=100;//xcf[ha_cgeshock[tindx1].ExoIndx];
               ha_cofvar1[tindx1].csolpupd=0;
               //if(rank==0&&i==0)printf("i %ld var %lf xsf %lf\n",tindx1,ha_cofvar1[tindx1].csolpupd,xcf[tindx1]);
            }
          }

          if(nohsl)MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time1);
          CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Last Update time %f\n",time1-time0);
          CHKERRQ(ierr);
          }

          
          if(!(subindx==subints-1&&sol==maxsol-1)) {
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxcf");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
              //return 1;
            }
            fwrite(xcf, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            free(xcf);
            xcf=NULL;
            //xcf=realloc (xcf,1*sizeof(ha_cgetype));
            
            if(subindx>0&&sol==0){            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxcO");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
              //return 1;
            }
            fwrite(xc0, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            }
            free(xc0);
            xc0=NULL;
            //xc0=realloc (xc0,2*sizeof(ha_cgetype));
//             if(subindx>0&&sol==0){            strcpy(tempfilenam,temdir);
//             strcat(tempfilenam,"_tempxcO12");
//             sprintf(tempchar, "%d",rank);
//             strcat(tempfilenam,tempchar);
//             strcat(tempfilenam,".bin");
//             if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
//               printf("Error opening file\n");
//               //return 1;
//             }
//             fwrite(xc012, sizeof(ha_cgetype),nvarele, tempvar);
//             fclose(tempvar);
//             }
//             xc012=realloc (xc012,2*sizeof(ha_cgetype));
// 
//             if(subindx>0&&sol==1){            strcpy(tempfilenam,temdir);
//             strcat(tempfilenam,"_tempxcO24");
//             sprintf(tempchar, "%d",rank);
//             strcat(tempfilenam,tempchar);
//             strcat(tempfilenam,".bin");
//             if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
//               printf("Error opening file\n");
//               //return 1;
//             }
//             fwrite(xc024, sizeof(ha_cgetype),nvarele, tempvar);
//             fclose(tempvar);
//             }
//             xc024=realloc (xc024,2*sizeof(ha_cgetype));
          }

    if(rank==rank_hsl&&sol==maxsol-1) {
            if(subindx==0){
              xc124=realloc (xc124,nvarele*sizeof(int));
              //memset(xc124,0,nvarele*sizeof(int));
              for(i=0; i<nvarele; i++)xc124[i]=6;
            }else{
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc124");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            xc124=realloc (xc124,nvarele*sizeof(int));
            freadresult=fread(xc124, sizeof(int),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            }

      for(i=0; i<nvarele; i++){
        j=0;
        if(xc12[i]>0){
        while (xc12[i] >= 10){
          xc12[i] /= 10;
          j++;
        }
        }else{
        while (xc12[i] <= -10){
          xc12[i] /= 10;
          j++;
        }
        }
        xc24[i]/=pow(10,j);
        j=abs(floor((xc12[i]-xc24[i])*100000));
        //if(rank==0)printf("12 %lf 24 %lf pres %ld\n",xc12[i],xc24[i],j);
        //if(j==0){
          //if(xc124[i]<6)xc124[i]=6;
        if(j!=0){//}else {
          if(j<10){
            if(xc124[i]>5)xc124[i]=5;
          }else{
            if(j<100){
              if(xc124[i]>4)xc124[i]=4;
            }else {
              if(j<1000){
                if(xc124[i]>3)xc124[i]=3;
              }else {
                if(j<10000){
                  if(xc124[i]>2)xc124[i]=2;
                }else{
                  if(xc124[i]>1)xc124[i]=1;
                }
              }
            }
          }
        }
      }
      if(subindx!=subints-1){
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc124");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
              //return 1;
            }
            fwrite(xc124, sizeof(int),nvarele, tempvar);
            fclose(tempvar);
            free(xc124);
            xc124=NULL;
            //xc124=realloc (xc124,1*sizeof(int));
      }
    }
            
          if(!(subindx==subints-1&&sol==maxsol-1)) {
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc12");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
              //return 1;
            }
            fwrite(xc12, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            free(xc12);
            xc12=NULL;
            //xc12=realloc (xc12,1*sizeof(ha_cgetype));
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc24");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
              //return 1;
            }
            fwrite(xc24, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            free(xc24);
            xc24=NULL;
            //xc24=realloc (xc24,1*sizeof(ha_cgetype));
          }
          
          
          strcpy(solchar,temdir);
          if(subindx<10)strcat(solchar,"_tempsol0");
          else strcat(solchar,"_tempsol");
          sprintf(tempchar, "%d", subindx);
          strcat(solchar,tempchar);
          sprintf(tempchar, "%d", sol);
          strcat(solchar,tempchar);
          strcat(solchar,".bin");
          printf("solchar %s\n",solchar);
          if ( (solution = fopen(solchar, "wb")) == NULL ) {
            printf("Error opening file\n");
            return 1;
          }
          fwrite(xcf, sizeof(ha_cgetype),nvarele, solution);
          fclose(solution);
        }
        free(x1);
        x1=NULL;
        //x1=realloc (x1,sizeof(ha_cgetype));
      }
    }
    long int *precis= (long int *) calloc (6,sizeof(long int));
    
    if(rank==rank_hsl) {
      for(i=0; i<nvarele; i++){
   switch(xc124[i]) {
      case 6:
         precis[5]+=1;
         break;
      case 5:
         precis[4]+=1;
         break;
      case 4:
         precis[3]+=1;
         break;
      case 3:
         precis[2]+=1;
         break;
      case 2:
         precis[1]+=1;
         break;
      default :
         precis[0]+=1;
   }
    }
    if(rank==0)printf("Accurate at 6 digits        %ld\nAccurate at 5 digits        %ld\nAccurate at 4 digits        %ld\nAccurate at 3 digits        %ld\nAccurate at 2 digits        %ld\nAccurate at 1 digit or none %ld\n",precis[5],precis[4],precis[3],precis[2],precis[1],precis[0]);
    }
    free(precis);
    xc0=realloc (xc0,sizeof(ha_cgetype));
    free(xc0);
    xc0=NULL;
    //xc012=realloc (xc012,sizeof(ha_cgetype));
    //xc024=realloc (xc024,sizeof(ha_cgetype));
    free(xc12);
    xc12=NULL;
    //xc12=realloc (xc12,sizeof(ha_cgetype));
    free(xc24);
    xc24=NULL;
    //xc24=realloc (xc24,sizeof(ha_cgetype));
    free(xc124);
    xc124=NULL;
    //xc124=realloc (xc124,sizeof(int));
    free(clag1);
    free(varchange);
    //time(&timeend);
    gettimeofday(&endtime, NULL);
    if(rank==0)printf("Mmid method calculation time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
    //if(rank==0)printf("Gragg method calculation time %f\n",difftime(timeend,timemulti));
              free(counteqs);
              free(counteqnoadds);
              free(countvarintra1s);
  }
    
  if(solmethod==20){
    MPI_Barrier(PETSC_COMM_WORLD);
  for(j=0;j<StoIter;j++){
    if(j>0){
              if(nohsl) {
                VecCreate(PETSC_COMM_WORLD,&vece);
              }
              else {
                VecCreate(PETSC_COMM_SELF,&vece);
              }
              if(nohsl) {
                VecSetType(vece,VECMPI);
              }
              else {
                VecSetType(vece,VECSEQ);
              }
              if(nesteddbbd==1)VecSetSizes(vece,localsize,VecSize);
              else VecSetSizes(vece,PETSC_DECIDE,VecSize);
              VecSetOption(vece, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
    }
    //Johansen(nohsl,VecSize,A,dnz,dnnz,onz,onnz,B,dnzB,dnnzB,onzB,onnzB,vecb,vece,rank,rank_hsl,mpisize,tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,&ha_cofvar,ncofele+nvarele,ncofele,nvarele,&ha_cgeshock,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,ha_eqadd,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,&xcf);
    //if(j==0)
      if(rank==rank_hsl)hnew_biupd(rank,tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,nvarele,10*laA,subints,1,1,0);
      
      ModMidPoint(nohsl,VecSize,&A,dnz,dnnz,onz,onnz,&B,dnzB,dnnzB,onzB,onnzB,&vecb,&vece,rank,rank_hsl,mpisize,tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,&ha_cofvar,ncofele+nvarele,ncofele,nvarele,&ha_cgeshock,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,ha_eqadd,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,subints,fcomm,&xcf,0);
    //else ModMidPoint(nohsl,VecSize,&A,dnz,dnnz,onz,onnz,&B,dnzB,dnnzB,onzB,onnzB,&vecb,&vece,rank,rank_hsl,mpisize,tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,&ha_cofvar,ncofele+nvarele,ncofele,nvarele,&ha_cgeshock,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,ha_eqadd,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,subints,fcomm,&xcf,0);
  printf("ncof %ld rank %d\n",ncof,rank);
//     if(rank==rank_hsl){
//       if(j==0)hnew_biupd(rank,tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,nvarele,10*laA,subints,1,1,0);
//       else hnew_biupd(rank,tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,nvarele,10*laA,subints,0,1,0);
//     }
  }
  }
    //if(rank==rank_hsl)hnew_biupd(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,laA);
  if(solmethod==21){
    MPI_Barrier(PETSC_COMM_WORLD);
    printf("Heere rank %d\n",rank);
      ModMidPoint(nohsl,VecSize,&A,dnz,dnnz,onz,onnz,&B,dnzB,dnnzB,onzB,onnzB,&vecb,&vece,rank,rank_hsl,mpisize,tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,&ha_cofvar,ncofele+nvarele,ncofele,nvarele,&ha_cgeshock,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,ha_eqadd,ndblock,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,begintime,subints,fcomm,&xcf,2);
  }
    
  if(rank==rank_hsl) {
    //for(i=0;i<VecSize;i++)printf("x0 %lf\n",xc[i]);
    //ofp = fopen("xvar1.txt", "w");
    //for(i=0; i<ncofele; i++) fprintf(ofp, "Cof Name: %s Cof Val %f\n",ha_cofvar[i].varname,xs[i]);
    //for(i=0; i<nvarele; i++) fprintf(ofp, "Var Name: %s Var Val %f Var change %f\n",ha_cofvar[ncofele+i].varname,xs[ncofele+i],xc[i]);
//for(i=0; i<nvar; i++) {
//  for(j=0; j<ha_var[i].matsize; j++) {
//    fprintf(ofp, "Var Name: %s Var Val %f Var change %f\n",ha_cofvar[ncofele+ha_var[i].begadd+j].varname,xs[ncofele+ha_var[i].begadd+j],xc[ha_cgeshock[ha_var[i].begadd+j].ExoIndx]);
//  }
//}
//for(i=0; i<nvarele; i++) {
//  if(ha_cgeshock[i].ShockId==0)
//    fprintf(ofp, "Var Name: %s Var Val %f Var change %f\n",ha_cofvar[ncofele+i].varname,xs[ncofele+i],xc[i]);
//  else
//    fprintf(ofp, "Var Name: %s Var Val %f Var change %f\n",ha_cofvar[ncofele+i].varname,xs[ncofele+i],xc[i]);
//}
//for(i=0; i<nvar; i++) printf("var %s change %d\n",ha_var[i].cofname,ha_var[i].begadd);
    //fclose(ofp);
    //printf("pc %lf vkb %lf\n",xc[106689],xc[106656]);
    for (i=niodata+noutdata; i<niodata+noutdata+nsoldata; i++) {
      //printf("logname %s filname %s\n",iodata[i].logname,iodata[i].filname);
      if (strcmp("solfiles",iodata[i].logname)==0) {
        strcpy(tempchar,iodata[i].filname);
        break;
      }
    }
    if(i==niodata+noutdata+nsoldata) {
      strcpy(tempchar,"solution");
    }
    strcpy(solchar,tempchar);
    strcat(solchar,".bin");
    printf("solchar %s\n",solchar);
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error opening file\n");
      return 1;
    }
    fwrite(xcf, sizeof(ha_cgetype),nvarele, solution);
    fclose(solution);
    //printf("i %ld xc %lf\n",nvarele-1,xc[nvarele-1]);
    strcpy(solchar,tempchar);
    strcat(solchar,".var");
    printf("solchar %s\n",solchar);
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error opening file\n");
      return 1;
    }
    fwrite(ha_var, sizeof(hcge_cof),nvar, solution);
    fclose(solution);
    strcpy(solchar,tempchar);
    strcat(solchar,".set");
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error opening file\n");
      return 1;
    }
    fwrite(ha_set, sizeof(ha_cgeset),nset, solution);
    fclose(solution);
    strcpy(solchar,tempchar);
    strcat(solchar,".sel");
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error opening file\n");
      return 1;
    }
    fwrite(ha_setele, sizeof(ha_cgesetele),nsetspace, solution);
    fclose(solution);
    uvadd modeldes[4];
    modeldes[0]=nsetspace;
    modeldes[1]=nvar;
    modeldes[2]=nvarele;
    modeldes[3]=(uvadd)nset;
    //printf("nset %d\n",nset);
    strcpy(solchar,tempchar);
    strcat(solchar,".mds");
    if ( (solution = fopen(solchar, "wb")) == NULL ) {
      printf("Error opening file\n");
      return 1;
    }
    fwrite(modeldes, sizeof(uvadd),4, solution);
    fclose(solution);
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  if(nowrites==0&&rank==0)for(i=0; i<noutdata; i++){
    hcge_wdata(tabfile,iodata[i+niodata].logname,iodata[i+niodata].filname,ha_set,nset,ha_setele,ha_cof,ncof,ncofele,ha_var,nvar,nvarele,ha_cofvar);
    printf("outfile %s\n",iodata[i+niodata].logname);
  }
  free(iodata);
  /*for(i=0; i<VecSize; i++) {
    VecSetValue(x,i,xc[i],INSERT_VALUES);
  }
  ierr = VecAssemblyBegin(x);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x);
  CHKERRQ(ierr);*/
  //PetscViewerASCIIOpen(PETSC_COMM_WORLD, "x2.txt", &viewer);
  //ierr = VecView(x,viewer);
  //CHKERRQ(ierr);
  //PetscViewerDestroy(&viewer);
  //printf("Hello world2 rank %d wr_id %d!\n",rank,wr_id);
  //if(rank==rank_hsl) if(remove("_tempvar.bin")==-1) {
  //printf("Unable to delete the file\n");
  //}
//hnew_update(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar1,ncofele+nvarele,ncofele);
//ofp = fopen("updated1.txt", "w");
//for(i=1; i<ncofele; i++) {
//  fprintf(ofp, "Cof Name: %s Cof Val %f\n",ha_cofvar[i].varname,ha_cofvar[i].varval);
//}
//fclose(ofp);
  free(countvarintra1);
  //free(ha_eqint);
  free(ha_eqadd);
  //free(ha_eqtime);
  //free(ha_eqreg);
  free(ndbbddrank1);
  free(counteq);
  free(counteqnoadd);
  free(ha_set);
  free(ha_setele);
  free(ha_cof);
  free(ha_var);
  free(ha_cgeshock);
  free(ha_cofvar);
  free(xcf);
//   if(solmethod==21){
//   char filename1[1024],j1name[1024];
//       if(rank<10)strcpy(j1name,"000");
//       if(rank<100&&rank>9)strcpy(j1name,"00");
//       if(rank<1000&&rank>99)strcpy(j1name,"0");
//       if(rank>=1000)j1name[0]='\0';
//       sprintf(filename1, "%d",rank);
//       strcat(j1name,filename1);
//       strcpy(filename1,"_biupd");
//       strcat(filename1,j1name);
//       strcat(filename1,".bin");
//       remove(filename1);
//   }
  //free(xs);
  //free(x1);
  printf("Hello world2! sof int %ld sof PetscInt %ld rank %d rankhsl %d\n",sizeof(int),sizeof(PetscInt),rank,rank_hsl);
  //free(b1);
  if(x0!=NULL)free(x0);
//**************************************************************************************
//**************************************END HSL*****************************************
//**************************************************************************************

  ierr = PetscFree(dnnzB);
  CHKERRQ(ierr);
  ierr = PetscFree(onnzB);
  CHKERRQ(ierr);
  ierr = PetscFree(dnnz);
  CHKERRQ(ierr);
  ierr = PetscFree(onnz);
  CHKERRQ(ierr);
  //ierr = VecDestroy(&b);
//CHKERRQ(ierr);
//ierr = VecDestroy(&b1);
  //CHKERRQ(ierr);
  //ierr = VecDestroy(&e);
  //CHKERRQ(ierr);
//ierr = VecDestroy(&e1);
//CHKERRQ(ierr);
  //ierr = VecDestroy(&x);
  //CHKERRQ(ierr);
  //ierr = MatDestroy(&A);
  //CHKERRQ(ierr);
  //ierr = MatDestroy(&B);
  //CHKERRQ(ierr);
  //ierr = MatDestroy(&G);
  //CHKERRQ(ierr);
  MPI_Comm_free(&HA_COMM);
  MPI_Comm_free(&HA1_COMM);
  ierr = PetscFinalize();
  CHKERRQ(ierr);
//**************************************************************************************
//****************************** END MATRIX FROM FORMULA********************************
//**************************************************************************************

//**************************************************************************************
//*********************************CLOSURE AND SHOCKS***********************************
//**************************************************************************************
  /*  char *closure="sj.cls";//, *shock="shock.shf";
    FILE * filehandle;
    char line[TABREADLINE],linecopy[TABREADLINE];
    int nexo;
    ha_cgeexovar *ha_exov= (ha_cgeexovar *) calloc ((nvarele-neq),sizeof(ha_cgeexovar));
    strcpy(commsyntax,"exogenous");
    filehandle = fopen(closure,"r");
    dcount=0;
    while (ha_cgertabl(commsyntax,filehandle,line))
    {
      while (ha_cgefrstr(line,"\n", ""));
      while (ha_cgedrcmt(line,"!"));
      while (ha_cgefrstr(line,"  ", " "));
      while (ha_cgefrstr(line," ;", ";"));
      while (ha_cgefrstr(line,"p_", ""));
      printf("line %s",line);
      strcpy(linecopy,line);
      nexo=ha_cgenfind(line, " ");
      //ha_cgevar *ha_exovar= (ha_cgevar *) calloc (nexo,sizeof(ha_cgevar));
      readitem = strtok(line," ");
      for (i=0;i<nexo;i++)
      {
        if (i==nexo-1) readitem = strtok(NULL,";");
        else readitem = strtok(line," ");
        //printf("readitem %s",readitem);
        while (ha_cgefrstr(readitem," ", ""));
          for (i=0;i<nvarele;i++)
          {
      j=0;
      while (readitem[j] != '\0') j++;
        //printf("readitem %s vaele %s\n",readitem,ha_cofvar[i+ncofele].varname);
            if (strncmp(ha_cofvar[i+ncofele].varname,readitem,j)==0)
            {
              strcpy(ha_exov[dcount].exoname,ha_cofvar[i+ncofele].varname);
              ha_exov[dcount].pos=i+ncofele;
              dcount++;
            }
          }
      }
    }
    if (dcount==nvarele-neq) printf("%s\n", "Closure is fine!");
    else printf("You need more %d %s\n", nvarele-neq-dcount,"exogenous variables!");
    //printf("dcout %d, neq %d, nvarele %d",dcount,neq,nvarele);
    for (i=0;i<dcount;i++) printf("exovar %s position %d\n",ha_exov[i].exoname,ha_exov[i].pos);
  */


//**************************************************************************************
//****************************** END CLOSURE AND SHOCKS********************************
//**************************************************************************************

  return 0;
}


