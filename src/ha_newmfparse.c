#include <ha_cgeglobal.h>

extern void spec51_rank_(int *INSIZE,ha_cgetype *cntl6,int *IRN, int *JCN, ha_cgetype *VA);
extern void spec51m_rank_(int *INSIZE,ha_cgetype *cntl6,int *IRN, int *JCN, ha_cgetype *VA,int *IRNA, int *JCNA, int *KEEP,ha_cgetype *w51, int *iw51);
//  extern void spec48_ssol_(forint *INSIZE,forint *IRN, forint *JCN, ha_cgetype *VA, ha_cgetype *B, ha_cgetype *X);
extern void spec48_ssol2la_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, ha_cgetype *B, ha_cgetype *X);
extern void spec48m_ssol2la_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, ha_cgetype *B, ha_cgetype *X);
extern void spec48_msol_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, ha_cgetype *B, ha_cgetype *X, PetscInt *IRNC, PetscInt *JCNC, PetscScalar *VAC,int *IRNB,int *JCNB,PetscScalar *VALUESB,ha_cgetype *VECBIVI,int *bivinzrow0,int *bivinzcol0);//, forint *IRNV, forint *JCNV, ha_cgetype *VAV
extern void spec48m_msol_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, ha_cgetype *B, ha_cgetype *X, PetscInt *IRNC, PetscInt *JCNC, PetscScalar *VAC,int *IRNB,int *JCNB,PetscScalar *VALUESB,ha_cgetype *VECBIVI,int *bivinzrow0,int *bivinzcol0);//, forint *IRNV, forint *JCNV, ha_cgetype *VAV
extern void spec48_esol_(int *INSIZE,int *IRN, ha_cgetype *VA,int *KEEP, ha_cgetype *B, ha_cgetype *X);
extern void spec48m_esol_(int *INSIZE,int *IRN, ha_cgetype *VA,int *KEEP, ha_cgetype *B, ha_cgetype *X);
extern void spec48_rpesol_(int *INSIZE,int *IRN, ha_cgetype *VA,int *KEEP, ha_cgetype *B, ha_cgetype *X,ha_cgetype *cntl,ha_cgetype *rinfo,ha_cgetype *error1,int *icntl,int *info,ha_cgetype *w,int *iw);
extern void spec48m_rpesol_(int *INSIZE,int *IRN, ha_cgetype *VA,int *KEEP, ha_cgetype *B, ha_cgetype *X,ha_cgetype *cntl,ha_cgetype *rinfo,ha_cgetype *error1,int *icntl,int *info,ha_cgetype *w,int *iw);
extern void spec48_single_(forint *indata,int *irn, int *jcn,ha_cgetype *b1, ha_cgetype *values,ha_cgetype *x1, int *neleperrow,int *ai1, MPI_Fint *fcomm);
extern void spec48_nomc66_(forint *indata, int *jcn,ha_cgetype *b1, ha_cgetype *values,ha_cgetype *x1, int *neleperrow, MPI_Fint *fcomm,forint *rowptrin, forint *colptrin);
extern void my_spar_add_(ha_cgetype *vecbivi, int *biviindx,int *nz1,ha_cgetype *vecbivi0,int *biviindx0,int *nz0,int *nz2);
extern void my_spar_addl_(ha_cgetype *vecbivi, long int *biviindx,int *nz1,ha_cgetype *vecbivi0,long int *biviindx0,int *nz0,int *nz2);
extern void my_spar_add3l_(ha_cgetype *vecbivi, long int *biviindx,long int *nz1,ha_cgetype *vecbivi0,long int *biviindx0,long int *nz0,long int *nz2);
extern void my_spar_add1_(ha_cgetype *vecbivi, int *biviindx,int *irn, int *jcn,int *nz1,ha_cgetype *vecbivi0,int *biviindx0,int *nz0,int *nz2,int *ncol);
extern void my_spar_add1l_(ha_cgetype *vecbivi, long int *biviindx,int *irn, int *jcn,int *nz1,ha_cgetype *vecbivi0,long int *biviindx0,int *nz0,int *nz2,int *ncol);
extern void my_spar_add4l_(ha_cgetype *vecbivi, long int *biviindx,int *irn, int *jcn,long int *nz1,ha_cgetype *vecbivi0,long int *biviindx0,long int *nz0,long int *nz2,int *ncol);
extern void my_spar_add2_(ha_cgetype *vecbivi, int *biviindx,int *irn, int *jcn,int *nz1,ha_cgetype *vecbivi0,int *biviindx0,int *nz0,int *nz2,int *ncol,ha_cgetype *vecbivi2,int *irn2, int *jcn2,int *j2,ha_cgetype *cntl3);
extern void my_spar_comp_(int *biviindx,int *nz1,int *biviindx0,int *nz0,int *nz2);
extern void my_spar_compl_(long int *biviindx,long int *nz1,long int *biviindx0,long int *nz0,long int *nz2);
//  extern void my_vec_comz_(ha_cgetype vecbivi,int *biviindx,int *col, int *row, int *colsize,int *nz0,int *nz1);
extern void prep48_alu_(int *INSIZE,int *IRN,int* JCN,ha_cgetype *VA);
extern void prep48_alu1_(int *INSIZE,int *IRN,int* JCN,ha_cgetype *VA,ha_cgetype *W,int *IW,int *KEEP);
extern void prep48_msol_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, PetscInt *IRNC, PetscInt *JCNC, PetscScalar *VAC,int *IRNB,int *JCNB,PetscScalar *VALUESB,ha_cgetype *VECBIVI,long int *bivinzrow0,int *bivinzcol0);//, forint *IRNV, forint *JCNV, ha_cgetype *VAV
extern void prep48m_msol_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, PetscInt *IRNC, PetscInt *JCNC, PetscScalar *VAC,int *IRNB,int *JCNB,PetscScalar *VALUESB,ha_cgetype *VECBIVI,long int *bivinzrow0,int *bivinzcol0,int* jcnb1,ha_cgetype *sol,ha_cgetype *b,ha_cgetype *w,int* iw,int* keep);
extern void spar_mulmin_(ha_cgetype* sol,int* nrow,int* nz,int* irn,int* jcn,ha_cgetype* va,ha_cgetype* res);
extern void spar_muladd_(ha_cgetype* sol,int* nrow,int* nz,int* irn,int* jcn,ha_cgetype* va,ha_cgetype* res);
extern void spar_mulnoadd_(ha_cgetype* sol,int* nrow,int* nz,int* irn,int* jcn,ha_cgetype* va,ha_cgetype* res);
extern void spar_vbiviadd_(ha_cgetype* sol,int* bvcol,long int* bvrow,long int* bvsize,int* nrow,int *ncol,int* nz,int* irn,int* jcn,ha_cgetype* va,ha_cgetype* res);
extern void patio_mat_(int* insizeda,int* IRN,int* JCN,ha_cgetype* VBIVI,int* IRN1A,int* JCN1A);

int HaNewMatVal(char *fname, char *commsyntax,ha_cgeset *ha_set,uvadd nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele,ha_cgeexovar *ha_cgeshock,uvadd ndblock,uvadd alltimeset,uvadd allregset,uvadd *ha_eqadd,uvadd *counteq,uvadd nintraeq,Mat A,Mat B) {
  FILE * filehandle;
  char tline[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],leftline[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[TABREADLINE],sumsyntax[NAMESIZE],lintmp[TABREADLINE];//,*p1=NULL;
  char *readitem=NULL,*p=NULL,*p1=NULL;//,*p2=NULL,*varpnts;
  ha_cgetype vval;
  PetscScalar    vvalpetsc;
  PetscInt Iindx,Jindx=0,Istart1,Iend1,matrow,rowindx;//,jfix=0
  PetscErrorCode ierr;
  ha_cgetype zerodivide=0;
  PetscMPIInt  mpisize1;
  bool isinproc;
  ierr = MatGetOwnershipRange(A,&Istart1,&Iend1);
  MPI_Comm_size(PETSC_COMM_WORLD,&mpisize1);
  //printf("Istart1 %d Iend1 %d\n",Istart1,Iend1);
  CHKERRQ(ierr);
  uvdim fdim,np,dcount,fdimlin=0,i4,sup,supset[MAXSUPSET];
  int totalsum,sumcount=1,sumcount1=0,lvar,lvar1,lvar2,lvar3,lvar4;
  uvadd lj,l1,i1=0,sumbegadd,dcountdim1[4*MAXVARDIM],dcountdim2[4*MAXVARDIM],dcountdim3[4*MAXVARDIM],nloops,nloopslin,nloopsfac,li3,nsumele,nsumele1,l2,eqindx=0,ltime;//,sizelinvars,totlinvars,templinvars
  int sumindx,npow,npar,nmul,nplu,ndiv,nmin,ha_calvarsize=0,nlinvars,leadlag,varindx1,varindx2;
  //uvadd *counteq1= (uvadd *) calloc (ndblock,sizeof(uvadd));
  uvadd j,sj,l,i3,i5,i,arsetdim=0,ha_calvardim=0;
  ha_cgesetindx *arSet1=NULL;
  hcge_calvars *ha_calvar1= NULL;
  //for(i=0; i<ndblock; i++)counteq1[i]=counteq[i];

  filehandle = fopen(fname,"r");
  matrow=0;//nintraeq;

  while (ha_cgertabl1(commsyntax,filehandle,line,ha_cofvar,ha_cof,ncof,&zerodivide,TABREADLINE)) {
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"line %s\n",line);
    if (strstr(line,"(default")==NULL) {
      ha_cgefrstr1(line, commsyntax, "");
      ha_cgefrstr1(line, "(linear)", "");
      while (ha_cgefrstr(line,"  ", " "));
      while (ha_cgefrchr(line, '[', '('));
      while (ha_cgefrchr(line, ']', ')'));
      while (ha_cgefrchr(line, '{', '('));
      while (ha_cgefrchr(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=ha_cgenfind(line, "(all,");
      //printf("row %d line %s\n",matrow,line);
      if (fdim==0) {
        readitem = strtok(line+1," ");
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      else {
        i=ha_cgerevfind(line, "(all,");
        readitem=line+i;
        readitem = strtok(readitem,")");
        readitem = strtok(NULL,"=");
        //printf("left %s\n",readitem);
        strcpy(vname,readitem);
        //printf("left %s linecopy %s\n",linecopy1,linecopy);
        strcpy(line,linecopy);
        //printf("left %s\n",line);
        readitem = strtok(line,"=");
        //printf("read %s\n",readitem);
        readitem = strtok(NULL,";");
        //printf("left %s\n",readitem);
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      while (ha_cgefrstr(readitem," ", ""));
      while (ha_cgerecovar(readitem)==1);
      hnew_intrpl(readitem);
      strcpy(tline,readitem);
      //ierr = PetscPrintf(PETSC_COMM_WORLD,"row %d: %s\n",matrow,readitem);
      strcpy(line1,readitem);
      npow=ha_cgenchf(readitem, '^');
      nmul=ha_cgenchf(readitem, '*');
      ndiv=ha_cgenchf(readitem, '/');
      nmul=nmul+ndiv;
      nplu=ha_cgenchf(readitem, '+');
      nmin=ha_cgenchf(readitem, '-');
      nplu=nplu+nmin;
      npar=ha_cgenchf(readitem, '(');

      strcpy(line,line1);
      readitem=line;
      //printf("read %s\n",readitem);
      np=ha_cgenfind(readitem,"p_");
      ha_calvardim=2*(npow+nmul+nplu+npar+1);
      hcge_calvars *ha_calvar= (hcge_calvars *) calloc (ha_calvardim,sizeof(hcge_calvars));
      hcge_linvars *LinVars= (hcge_linvars *) calloc (np,sizeof(hcge_linvars));
      //unsigned long int *LinVarsBA= (unsigned long int *) calloc (np,sizeof(unsigned long int));
      i3=0;
      //sizelinvars=0;
      lvar=0;
      //printf("linecopy %s\n",linecopy);
      for (i=0; i<np; i++) {
        varindx2=0;
        while(-1<0) {
          varindx1=ha_cgefind(readitem+varindx2,"p_");
          if(varindx1==-1) break;
          varindx2=varindx2+varindx1;
          if(varindx2==0||readitem[varindx2-1]=='*'||readitem[varindx2-1]=='+'||readitem[varindx2-1]=='-'||readitem[varindx2-1]=='('||readitem[varindx2-1]==',') break;
          else varindx2++;
        }
        if(varindx1==-1) break;
        lvar=lvar+varindx2;
        readitem=readitem+varindx2+2;

        //lvar=lvar+ha_cgefind(readitem,"p_");
        //readitem=readitem+ha_cgefind(readitem,"p_")+2;
        p=strpbrk(readitem,"{+*-/^)");
        if(p!=NULL) {
          strncpy(vname,readitem,p-readitem);
          vname[p-readitem]='\0';
          //printf("p %s np %d\n",vname,np);
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(ha_var[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          p=strpbrk(readitem,"}+*-/^)");
          if (*p=='}') {
            p++;
          }
          strncpy(tline,readitem,p-readitem);
          tline[p-readitem]='\0';
          switch (ha_var[l].size) {
          case 0:
            break;
          case 1:
            p = strtok(tline,"{");
            p = strtok(NULL,"}");
            leadlag=0;
            hnew_arset(p,&leadlag);
            strcpy(LinVars[i3].dimnames[0],p);
            LinVars[i3].dimleadlag[0]=leadlag;
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=ha_cgefind(linecopy,lintmp);
            //printf("lin %s\n",lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[0],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              //printf("lin1 %s\n",lintmp);
              lvar1=ha_cgenfind(linecopy,lintmp);
              lvar3=ha_cgefind(linecopy,lintmp);
              //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
                  if (lvar4>-1&&lvar4<lvar) {
                    lvar3=lvar3+lvar4+4;
                  }
                  else {
                    break;
                  }
                }
              //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
              p1=&linecopy[0]+lvar3;
              //p1=&linecopy[0]+ha_cgefind(linecopy,lintmp);
              p1=p1+strlen(lintmp);
              //printf("lin %s\n",lintmp);
              strncpy(LinVars[i3].dimsetnames[0],p1,strchr(p1,',')-p1);
            }
            //printf("vname %s lin %s\n",vname,LinVars[i3].dimsetnames[0]);
            //printf("lin %s indx %s\n",ha_var[l].cofname,p);
            break;
          default:
            p = strtok(tline,"{");
            for (i4=0; i4<ha_var[l].size-1; i4++) {
              p = strtok(NULL,",");
              leadlag=0;
              hnew_arset(p,&leadlag);
              strcpy(LinVars[i3].dimnames[i4],p);
              LinVars[i3].dimleadlag[i4]=leadlag;
              strcpy(lintmp,"(all,");
              strcat(lintmp,p);
              strcat(lintmp,",");
              l1=ha_cgefind(linecopy,lintmp);
              //printf("lin %s\n",lintmp);
              if (l1>-1) {
                p1=&linecopy[0]+l1;
                strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
              }
              else {
                strcpy(lintmp,"sum(");
                strcat(lintmp,p);
                strcat(lintmp,",");
                //printf("lin1 %s\n",lintmp);
                lvar1=ha_cgenfind(linecopy,lintmp);
                lvar3=ha_cgefind(linecopy,lintmp);
                //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
                if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                    lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
                    if (lvar4>-1&&lvar4<lvar) {
                      lvar3=lvar3+lvar4+4;
                    }
                    else {
                      break;
                    }
                  }
                //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
                p1=&linecopy[0]+lvar3;
                //p1=&linecopy[0]+ha_cgefind(linecopy,lintmp);
                p1=p1+strlen(lintmp);
                //printf("lin %s p1 %s len %d\n",lintmp,p1,strchr(p1,',')-p1-strlen(lintmp));
                strncpy(LinVars[i3].dimsetnames[i4],p1,strchr(p1,',')-p1);
                //printf("lin %s\n",LinVars[i3].dimsetnames[i4]);
              }
              //printf("vname %s lin %s\n",vname,LinVars[i3].dimsetnames[i4]);
            }
            p = strtok(NULL,"}");
            leadlag=0;
            hnew_arset(p,&leadlag);
            strcpy(LinVars[i3].dimnames[i4],p);
            LinVars[i3].dimleadlag[i4]=leadlag;
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=ha_cgefind(linecopy,lintmp);
            //printf("lin %s\n",lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=ha_cgenfind(linecopy,lintmp);
              lvar3=ha_cgefind(linecopy,lintmp);
              //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
                  if (lvar4>-1&&lvar4<lvar) {
                    lvar3=lvar3+lvar4+4;
                  }
                  else {
                    break;
                  }
                }
              //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
              p1=&linecopy[0]+lvar3;
              p1=p1+strlen(lintmp);
              //printf("lin %s p1 %s\n",lintmp,p1);
              strncpy(LinVars[i3].dimsetnames[i4],p1,strchr(p1,',')-p1);
            }
            //printf("vname %s lin %s l1 %d\n",vname,LinVars[i3].dimsetnames[i4],l1);
            break;
          }
          i3++;
          //}
          //LinVars[i][p-readitem]='\0';
        }
        else {
          strcpy(vname,readitem);
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(ha_var[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          i3++;
          //}
          //break;
        }
      }
      nlinvars=i3;
      i3=0;
      for (i=0; i<nlinvars; i++) i3+=ha_var[LinVars[i].LinVarIndx].size;
      arsetdim=fdim+i3;
      ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (arsetdim,sizeof(ha_cgesetindx));

      strcpy(line,linecopy);
      //printf("adim1 %ld fdim1 %d\n",arsetdim,fdim);
      nloops=1;
      if (fdim>0) {
        for (i=0; i<fdim; i++) {
          if(i==0) {
            readitem = strtok(line,",");
          }
          else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
          strcpy(arSet[i].arIndx,readitem);
          readitem = strtok(NULL,")");
          //strcpy(arSet[i].arSet,readitem);
          for (i4=0; i4<nset; i4++) if(strcmp(readitem,ha_set[i4].setname)==0) {
              arSet[i].setid=i4;
              //arSet[i].SetSize=ha_set[i4].size;
              //arSet[i].SetBegAdd=ha_set[i4].begadd;
              //arSet[i].subsetid=ha_set[i4].subsetid;
              //if(ha_set[i4].subsetid==1) {
              //  arSet[i].SuperSetSize=ha_set[i4].supersetsize;
              //check=true;
              //} else {
              //  arSet[i].SuperSetSize=ha_set[i4].size;
              //}
              break;
            }
          nloops=nloops*ha_set[arSet[i].setid].size;
        }
        dcountdim1[fdim-1]=1;
        for (dcount=fdim-2; dcount>-1; dcount--) {
          dcountdim1[dcount]=ha_set[arSet[dcount+1].setid].size*dcountdim1[dcount+1];
        }
        /*for (dcount=ha_eqtime[eqindx]-1; dcount>-1; dcount--) {
          dcountdim3[dcount]=dcountdim1[dcount]/arSet[ha_eqtime[eqindx]].SetSize;
        }
        for (dcount=fdim-1; dcount>ha_eqtime[eqindx]-1; dcount--) {
          dcountdim3[dcount]=dcountdim1[dcount];
        }*/
      }
      isinproc=false;
      if(mpisize1>1) {
        for (lj=0; lj<nloops; lj++) {
          Jindx=ha_eqadd[matrow+lj];
          if(Jindx>=Istart1&&Jindx<Iend1) {
            isinproc=true;
            break;
          }
        }
      }
      else isinproc=true;
      //if(eqindx==0)for(i=0;i<fdim;i++)printf("eqtime %d dc3 %d dc1 %d\n",ha_eqtime[eqindx],dcountdim3[i],dcountdim1[i]);
      if(isinproc) {
        strcpy(line,line1);
        readitem=line;
        strcpy(sumsyntax,"sum(");
        totalsum=hcge_nsum(readitem,sumsyntax);
        hcge_sumcof *sum_cof= (hcge_sumcof *) calloc (totalsum*nlinvars+1,sizeof(hcge_sumcof));
        sumcount=0;
        while (hcge_dsum(readitem,sumsyntax,sum_cof,arSet,ha_set,nset,fdim+1,sumcount)==1) {
          sumcount++;
        }
        totalsum=sumcount;
        //for (i3=0; i3<totalsum; i3++) printf("varsum %d %s\n",i3,sum_cof[i3].sumname);
        li3=0;
        for (i4=0; i4<totalsum; i4++) {
          i1=1;
          for(sj=0; sj<sum_cof[i4].size; sj++) {
            i1=i1*ha_set[sum_cof[i4].setid[sj]].size;
          }
          sum_cof[i4].begadd=li3;
          li3=li3+i1;
        }
        nsumele=li3;
        //ha_cgemsumele *outlin_sumele= (ha_cgemsumele *) calloc (olinnsumele+1,sizeof(ha_cgemsumele));
        ha_cgesumele *ha_sumele= (ha_cgesumele *) calloc (nsumele*nlinvars,sizeof(ha_cgesumele));

        //printf("OK!!! sumele %d vars %d e %f\n",nsumele,nlinvars,ha_sumele[0].varval);
        strcpy(line,line1);
        readitem=line;
        sumcount=0;
        while (hlin_dsum(readitem,sumsyntax,sum_cof,arSet,ha_set,nset,fdim+1,sumcount)==1) {
          sumcount++;
        }
        //printf("OK1!!!\n");
        totalsum=sumcount;
        //for (i3=0; i3<totalsum; i3++) printf("sumvar %s\n",sum_cof[i3].sumname);
        li3=0;
        for (i4=0; i4<totalsum; i4++) {
          i1=1;
          for(sj=0; sj<sum_cof[i4].size; sj++) {
            i1=i1*ha_set[sum_cof[i4].setid[sj]].size;
          }
          sum_cof[i4].begadd=li3;
          li3=li3+i1;
        }
        nsumele=li3;
        for (i3=0; i3<totalsum; i3++) {
          i1=1;
          sum_cof[i3].antidims[sum_cof[i3].size-1]=1;
          for(sj=sum_cof[i3].size-2; sj>-1; sj--) {
            sum_cof[i3].antidims[sj]=sum_cof[i3].antidims[sj+1]*ha_set[sum_cof[i3].setid[sj+1]].size;
          }
        }
        //printf("OK2!!!\n");
        //hcge_rsumele(sum_cof,totalsum,ha_set,nset,ha_setele,nsumele);
        //if(totalsum>0){
        sumcount=0;
        strcpy(line,line1);
        readitem=line;
        sumindx=0;
        //printf("OK3!!!\n");
        while (hlin_calsum(readitem,sumsyntax,ha_set,nset,ha_setele,ha_cofvar,ncofvar,ncofele,ha_cof,ncof,ha_var,nvar,sum_cof,totalsum,ha_sumele,nsumele,ha_calvar,arSet,fdim+1,&sumindx,sumcount,zerodivide)==1) {
          sumcount++;
        }
        strcpy(line1,readitem);
        //}
        //printf("aline1 %s\n",line1);
        //for (i=0;i<nsumele;i++) printf("varele %s varval%f\n",ha_sumele[i].varname,ha_sumele[i].varval);
        sumbegadd=nsumele;
        //sumcount=0;
        sumcount1=sumcount;
        for (i=0; i<nlinvars; i++) {
          Jindx=ha_eqadd[matrow];
          if(Jindx>=Iend1)continue;
          i3=0;
          nloopslin=nloops;
          if (fdim==0) {
            for (l2=0; l2<ha_var[LinVars[i].LinVarIndx].size; l2++) {
              strcpy(arSet[i3].arIndx,LinVars[i].dimnames[l2]);
              //strcpy(arSet[i3].arSet,LinVars[i].dimsetnames[l2]);
              for (i4=0; i4<nset; i4++) if(strcmp(ha_set[i4].setname,LinVars[i].dimsetnames[l2])==0) {
                  break;
                }
              arSet[i3].setid=i4;
              //arSet[i3].SetSize=ha_set[i4].size;
              //arSet[i3].SetBegAdd=ha_set[i4].begadd;
              nloopslin=nloopslin*ha_set[arSet[i3].setid].size;
              //printf("set %s setsize %d nloop %d nlooplin %d\n",ha_set[i4].setname,arSet[i3].SetSize,nloops,nloopslin);
              i3++;
            }
          }
          else for (l2=0; l2<ha_var[LinVars[i].LinVarIndx].size; l2++) {
              for (sj=0; sj<fdim; sj++) {
                if(strcmp(arSet[sj].arIndx,LinVars[i].dimnames[l2])==0) {
                  break;
                }
                else if (sj==fdim-1) {
                  strcpy(arSet[fdim+i3].arIndx,LinVars[i].dimnames[l2]);
                  //strcpy(arSet[fdim+i3].arSet,LinVars[i].dimsetnames[l2]);
                  for (i4=0; i4<nset; i4++) if(strcmp(ha_set[i4].setname,LinVars[i].dimsetnames[l2])==0) {
                      break;
                    }
                  arSet[fdim+i3].setid=i4;
                  //arSet[fdim+i3].SetSize=ha_set[i4].size;
                  //arSet[fdim+i3].SetBegAdd=ha_set[i4].begadd;
                  nloopslin=nloopslin*ha_set[arSet[fdim+i3].setid].size;
                  i3++;
                }
              }
            }
          nloopsfac=(uvadd)nloopslin/nloops;
          //Jindx=ha_eqadd[matrow+(uvadd)lj/nloopsfac];
          //if(Jindx>=Istart1&&Jindx<Iend1)
          Jindx=ha_eqadd[matrow+(uvadd)(nloopslin-1)/nloopsfac];
          if(Jindx<Istart1)continue;
          
          fdimlin=fdim+i3;
          //printf("fdim %d fdimlin %d nloop %ld nlooplin %ld\n",fdim,fdimlin,nloops,nloopslin);
          //if (matrow==31680)
          //for (i4=0;i4<fdimlin;i4++) printf("indx1 %s set %s\n",arSet[i4].arIndx,ha_set[arSet[i4].setid].setname);
          if (i3>0) {
            dcountdim2[fdimlin-1]=1;
            //printf("set %s\n",arSet[fdimlin-1].arSet);
            for (dcount=fdimlin-2; dcount>-1; dcount--) {
              //printf("set %s\n",arSet[dcount].arSet);
              dcountdim2[dcount]=ha_set[arSet[dcount+1].setid].size*dcountdim2[dcount+1];
            }
          }
          else for (dcount=0; dcount<fdim; dcount++) {
              dcountdim2[dcount]=dcountdim1[dcount];
            }

          //int *LinVarsele= (int *) calloc (ha_var[LinVars[i].LinVarIndx].matsize,sizeof(int));
          //printf("aline1 %s\n",line1);
          strcpy(leftline,line1);
          hcge_rlinzero(leftline,i);
          //printf("linreplsum %s\n",leftline);
          hcge_repllin(leftline,0);
          //printf("linreplsum %s\n",leftline);
          strcpy(sumsyntax,"sum(");
          hlin_replsum(leftline,sumsyntax,i,LinVars,ha_var);
          //hcge_repllin(leftline,1);
          //printf("linreplsum1 %s\n",leftline);
          hlin_rlinone(leftline,LinVars,i,ha_var);
          //printf("line one %s\n",leftline);
          strcpy(line,leftline);
          readitem=line;
          //totalsum=hcge_nsum(readitem,sumsyntax);
          //hcge_sumcof *sum_cof= (hcge_sumcof *) calloc (totalsum+1,sizeof(hcge_sumcof));
          //printf("reada3 %s\n",readitem);
          sumcount1=sumcount;
          while (hcge_dsum(readitem,sumsyntax,sum_cof,arSet,ha_set,nset,fdimlin+1,sumcount1)==1) {
            sumcount1++;
          }
          //for (i4=0; i4<totalsum; i4++)for (i3=0; i3<sum_cof[i4].size; i3++) printf("1sumname %s sum set %s\n",sum_cof[i4].sumname,ha_set[sum_cof[i4].setid[i3]].setname);
          //printf("reada3 %s\n",readitem);
          totalsum=sumcount1;
          //if(i==0)i3=nsumele;
          //printf("sumcount %d total sum %d\n",sumcount,totalsum);
          for (i4=sumcount; i4<totalsum; i4++) {
            i1=1;
            for(sj=0; sj<sum_cof[i4].size; sj++) {
              i1=i1*ha_set[sum_cof[i4].setid[sj]].size;
              //printf("sum %s set %s\n",sum_cof[i4].sumname,ha_set[sum_cof[i4].setid[sj]].setname);
            }
            sum_cof[i4].begadd=sumbegadd;
            sumbegadd=sumbegadd+i1;
          }
          nsumele1=sumbegadd+nsumele;
          //printf("nsumele %d\n",nsumele);
          for (i3=sumcount; i3<totalsum; i3++) {
            i1=1;
            sum_cof[i3].antidims[sum_cof[i3].size-1]=1;
            for(sj=sum_cof[i3].size-2; sj>-1; sj--) {
              sum_cof[i3].antidims[sj]=sum_cof[i3].antidims[sj+1]*ha_set[sum_cof[i3].setid[sj+1]].size;
              //printf("sumanti %d\n",sum_cof[i].antidims[j]);
            }
          }
          //printf("totalsum %d\n",totalsum);
          //ha_cgemsumele *ha_sumele= (ha_cgemsumele *) calloc (nsumele+1,sizeof(ha_cgemsumele));
          //printf("totalsum %d\n",totalsum);
          //hcge_rsumele(sum_cof,totalsum,ha_set,nset,ha_setele,ha_sumele);
          //if (matrow==382) for(i3=0;i3<nsumele;i3++)printf("sumele %s\n",ha_sumele[i3].varname);
          //if (matrow==31680) printf("read %s\n",readitem);
          //sumcount1=sumcount;
          strcpy(line,leftline);
          readitem=line;
          //sumindx1=sumindx;
          //printf("i %d sumcount1 %d sumindx1 %d\n",i,sumcount1,sumindx1);
          //printf("reada %s\n",readitem);
          //for (i4=0; i4<totalsum; i4++)for (i3=0; i3<sum_cof[i4].size; i3++) printf("sumname %s sum set %s\n",sum_cof[i4].sumname,ha_set[sum_cof[i4].setid[i3]].setname);
          //if(totalsum>sumcount){
          //printf("read %s\n",readitem);
          while (hnew_calsum(readitem,sumsyntax,ha_set,nset,ha_setele,ha_cofvar,ncofvar,ncofele,ha_cof,ncof,ha_var,nvar,sum_cof,totalsum,ha_sumele,nsumele1,ha_calvar,arSet,fdimlin+1,&sumindx,sumcount,zerodivide)==1) {
            sumcount++;
          }
          //}
          //strcpy(line1,readitem);
          //printf("read1 %s\n",readitem);
          //if (matrow==32)for (j=0;j<sumindx;j++) printf("varele %s varval%d\n",sum_cof[j].sumname,sum_cof[j].begadd);
          //if (matrow==32) for (j=0;j<nsumele1;j++) printf("varval%f\n",ha_sumele[j].varval);
          //if(readitem[0]=='+'||readitem[0]=='-') {memmove(readitem+1,readitem,strlen(readitem)+1);readitem[0]='0';}
          //printf("read %s\n",readitem);
          //ierr = PetscPrintf(PETSC_COMM_WORLD,"read %s\n",readitem);
          ha_calvarsize=0;
          ha_newfparse(readitem,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,&ha_calvarsize,arSet,fdimlin);
          //for (lj=0; lj<ha_calvarsize; lj++) ierr = PetscPrintf(PETSC_COMM_WORLD,"j %d oper %d tmpvar %s tmpval var1 %d var2 %d\n",lj,ha_calvar[lj].Oper,ha_calvar[lj].TmpVarName,ha_calvar[lj].Var1Type,ha_calvar[lj].Var2Type);
          //printf("cal var size %d total size %d OK!!!\n",ha_calvarsize,2*(npow+nmul+nplu+npar+1));
          //ierr = PetscPrintf(PETSC_COMM_WORLD,"row %d read %s\n",Jindx,readitem);
          //printf("nloop %d\n",nloopslin);
          for (dcount=0; dcount<ha_var[LinVars[i].LinVarIndx].size; dcount++) {
            for (i4=0; i4<fdimlin; i4++) {
              if (strcmp(LinVars[i].dimnames[dcount],arSet[i4].arIndx)==0) {
                dcountdim3[dcount]=i4;
                break;
              }
            }
          }
          for(dcount=0; dcount<MAXSUPSET; dcount++)supset[dcount]=0;
          for (dcount=0; dcount<ha_var[LinVars[i].LinVarIndx].size; dcount++) {
            if(ha_set[ha_var[LinVars[i].LinVarIndx].setid[dcount]].size!=ha_set[arSet[dcountdim3[dcount]].setid].size) {
              //printf("nrow %d vaset %s set ste %d\n",matrow,ha_set[ha_var[LinVars[i].LinVarIndx].setid[dcount]].setname,ha_set[arSet[dcountdim3[dcount]].setid].subsetid[1]);
              for(sup=1; sup<MAXSUPSET; sup++)if(ha_var[LinVars[i].LinVarIndx].setid[dcount]==ha_set[arSet[dcountdim3[dcount]].setid].subsetid[sup]) {
                  supset[dcount]=sup;
                  break;
                }
            }
            else supset[dcount]=0;
          }
        #pragma omp parallel private(lj,Jindx,i3,sj,i5,l2,dcount,l1,li3,Iindx,ierr,arSet1,ha_calvar1,vval) shared(ha_cofvar,arSet)
        {
        //j8=0;
        if(omp_get_thread_num()!=0){
          arSet1=realloc(arSet1,arsetdim*sizeof(ha_cgesetindx));
          memcpy(arSet1,arSet,arsetdim*sizeof(ha_cgesetindx));
          ha_calvar1=realloc(ha_calvar1,ha_calvardim*sizeof(hcge_calvars));
          memcpy(ha_calvar1,ha_calvar,ha_calvardim*sizeof(hcge_calvars));
        }else{
          ha_calvar1=ha_calvar;
          arSet1=arSet;
        }
          ha_cgetype *value= (ha_cgetype *) calloc (nloopsfac,sizeof(ha_cgetype));
          PetscInt *jcn= (PetscInt *) calloc (nloopsfac,sizeof(PetscInt));
          ha_cgetype *valueb= (ha_cgetype *) calloc (nloopsfac,sizeof(ha_cgetype));
          PetscInt *jcnb= (PetscInt *) calloc (nloopsfac,sizeof(PetscInt));
          //printf("HHHH\n");
        #pragma omp for
          //for (lj=0; lj<nloopslin; lj++) {
          for (i5=0; i5<nloops; i5++) {
            Jindx=ha_eqadd[matrow+i5];//Jindx=ha_eqadd[matrow+(uvadd)lj/nloopsfac];
            if(Jindx>=Istart1&&Jindx<Iend1) {
              i3=0;
              sj=0;
              for (lj=i5*nloopsfac; lj<(i5+1)*nloopsfac; lj++) {
              //if(Jindx>92734)ierr = PetscPrintf(PETSC_COMM_WORLD,"read end %d\n",lj);
              //printf("lj %ld fdimlin %ld\n",lj,fdimlin);
              l2=lj;
              //rowindx=0;
              for (dcount=0; dcount<fdimlin; dcount++) {
                l1=(uvadd) l2/dcountdim2[dcount];
                arSet1[dcount].indx=l1;
                //if(ha_eqint[eqindx]) {
                //  if (dcount<fdim&&dcount!=ha_eqtime[eqindx]) {
                //    rowindx=rowindx+l1*dcountdim3[dcount];
                //  }
                //  if(dcount==ha_eqtime[eqindx])ltime=l1;
                //} else {
                //if (dcount<fdim)rowindx=rowindx+l1*dcountdim1[dcount];
                //}
                l2=l2-l1*dcountdim2[dcount];
              }
              li3=0;
              //printf("OKKKK\n");
              for (dcount=0; dcount<ha_var[LinVars[i].LinVarIndx].size; dcount++) {
                //printf("d %d l3 %d\n",dcount,(arSet[dcountdim3[dcount]].indx+LinVars[i].dimleadlag[dcount])*ha_var[LinVars[i].LinVarIndx].antidims[dcount]);
                //printf("d1 %d\n",supset[dcount]);
                //for (i4=0; i4<fdimlin; i4++) {
                //if (strcmp(LinVars[i].dimnames[dcount],arSet[i4].arIndx)==0) {
                if(supset[dcount]==0) {
                  li3=li3+(arSet1[dcountdim3[dcount]].indx+LinVars[i].dimleadlag[dcount])*ha_var[LinVars[i].LinVarIndx].antidims[dcount];
                }
                else {
                  //for(sup=1;sup<MAXSUPSET;sup++)if(ha_var[LinVars[i].LinVarIndx].setid[dcount]==ha_set[arSet[dcountdim3[dcount]].setid].subsetid[sup]){
                  li3=li3+(ha_setele[ha_set[arSet1[dcountdim3[dcount]].setid].begadd+arSet1[dcountdim3[dcount]].indx].setsh[supset[dcount]]+LinVars[i].dimleadlag[dcount])*ha_var[LinVars[i].LinVarIndx].antidims[dcount];
                  //break;
                  //}
                }
                //break;
                //}
                //}
              }
              //printf("OKKKK\n");
              //if(matrow==31680)ierr = PetscSynchronizedPrintf(PETSC_COMM_SELF,"read 1 %s\n",readitem);
//              if(Jindx==16504911){
//              vval=ha_newfpcalshow(ha_cofvar,ha_set,ha_setele,ha_sumele,ha_calvar,ha_calvarsize,arSet,fdimlin,zerodivide,1);
//                printf("Jindx1 %d Indx1 %d vval %lf cals %d\n",Jindx,Iindx,vval,ha_calvarsize);
//                for(dcount=0;dcount<ha_calvarsize;dcount++){
//                  printf("var 1 %d var1val %lf var2 %d var2val %lf\n",ha_calvar[dcount].Var1BegAdd,ha_calvar[dcount].Var1Val,ha_calvar[dcount].Var2BegAdd,ha_calvar[dcount].Var2Val);
//                }
//              } else
              vval=ha_newfpcal(ha_cofvar,ha_set,ha_setele,ha_sumele,ha_calvar1,ha_calvarsize,arSet1,fdimlin,zerodivide);
              //if(Jindx>79127)ierr = PetscSynchronizedPrintf(PETSC_COMM_SELF,"read 2 %s vval %f\n",readitem,vval);
              //if(ha_eqint[eqindx]) Jindx=counteq1[ha_setele[arSet[ha_eqtime[eqindx]].SetBegAdd+ltime].setsh]+rowindx;
              //else Jindx=matrow+rowindx;
              //printf("eq %ld\n",matrow+(uvadd)lj/nloopsfac);
              //if(Jindx==21890||Jindx==21891)ierr = PetscSynchronizedPrintf(PETSC_COMM_SELF,"read 2a %s vval %f\n",line1,vval);
              //if(Jindx==21891)for(i=0;i<fdim;i++)printf("eqtime %d dc3 %d dc1 %d\n",ha_eqtime[eqindx],dcountdim3[i],dcountdim1[i]);
              //printf("var %d\n",ha_var[LinVars[i].LinVarIndx].begadd+li3);
              Iindx=ha_cgeshock[ha_var[LinVars[i].LinVarIndx].begadd+li3].ExoIndx;
              //if(Jindx<432&&Iindx>432&&Iindx<1296)printf("line %s\n Linvar %s Jindx %d Iindx %d vval %lf\n",linecopy,ha_var[LinVars[i].LinVarIndx].cofname,Jindx,Iindx,vval);
              //printf("leq %d Jindx %d Iindx %d vval %lf\n",matrow+(uvadd)lj/nloopsfac,Jindx,Iindx,vval);
              //printf("lin %s l3 %d Iindx %d begadd %d\n",ha_var[LinVars[i].LinVarIndx].cofname,li3,Iindx,ha_var[LinVars[i].LinVarIndx].begadd);
              //if(Jindx==178)printf("row %d line %s\n",matrow,linecopy);
              //if(Jindx==26422667)printf("ltime %d rowindx %d Jindx %d lj %d timeindx %d fdimlin %d set %d vval %f\n",ltime,rowindx,Jindx,lj,ha_eqtime[eqindx],fdimlin,ha_setele[arSet[ha_eqtime[eqindx]].SetBegAdd+ltime].setsh,vval);
              //if(Jindx==26422667)for (dcount=0; dcount<fdim; dcount++)printf("dcount %d\n",dcountdim3[dcount]);
              //if(Jindx==31680) for (dcount=0; dcount<fdimlin; dcount++) printf("ar %s arindx %d vval %f\n",arSet[dcount].arIndx,arSet[dcount].indx,vval);
              //if(Jindx>92734)ierr = PetscPrintf(PETSC_COMM_WORLD,"read 2b %s vval %f\n",readitem,vval);
              //if(eqindx==0)printf("J %d I %d Istart %d Iend %d vval %f\n",Jindx,Iindx,Istart1,Iend1,vval);
              //if(Istart1<=Jindx&&Jindx<Iend1) {
              //if(Jindx==16504911)printf("Jindx1 %d Indx1 %d vval %lf\n",Jindx,Iindx,vval);
              //#pragma omp critical
              //{
              if (!ha_cgeshock[ha_var[LinVars[i].LinVarIndx].begadd+li3].ShockId&&vval!=0) {
                //vvalpetsc=(PetscScalar)vval;
                //#pragma omp atomic
                value[i3]=vval;
                //irn[i3]=Jindx;
                jcn[i3]=Iindx;
                i3++;
                //ierr = MatSetValue(A,Jindx,Iindx,(PetscScalar)vval,ADD_VALUES);//MatSetValues(A,1,&Jindx,1,&Iindx,&vvalpetsc,ADD_VALUES);
                //if(Jindx>=5819901&&Jindx<=5819905)printf("Jindx %d Indx %d vval %lf\n",Jindx,Iindx,vval);
                //CHKERRQ(ierr);
              }
              if (ha_cgeshock[ha_var[LinVars[i].LinVarIndx].begadd+li3].ShockId&&vval!=0) {
                //vval=-vval;
                //vvalpetsc=(PetscScalar)vval;
                //#pragma omp atomic
                valueb[sj]=-vval;
                jcnb[sj]=Iindx;
                sj++;
                //ierr = MatSetValue(B,Jindx,Iindx,(PetscScalar)vval,ADD_VALUES);//MatSetValues(B,1,&Jindx,1,&Iindx,&vvalpetsc,ADD_VALUES);
                //CHKERRQ(ierr);
                //if(Jindx>=16514276-1000&&Jindx<16514276)printf("Jindx %d Indx %d vval %lf\n",Jindx,Iindx,vval);
              }
              //}
              //if(Jindx>92734)ierr = PetscPrintf(PETSC_COMM_WORLD,"read 2c %s vval %f\n",readitem,vval);
            }
            if(i3>0)MatSetValues(A,1,&Jindx,i3,jcn,value,ADD_VALUES);
            if(sj>0)MatSetValues(B,1,&Jindx,sj,jcnb,valueb,ADD_VALUES);
            }
            //if(Jindx>92734)ierr = PetscPrintf(PETSC_COMM_WORLD,"read 3 %s\n",readitem);
          }
          free(value);
          free(jcn);
          free(valueb);
          free(jcnb);
        if(omp_get_thread_num()!=0){
          free(arSet1);
          arSet1=NULL;
          free(ha_calvar1);
          ha_calvar1=NULL;
        }else{
          ha_calvar1=NULL;
          arSet1=NULL;
        }
        }
          //ierr = PetscPrintf(PETSC_COMM_WORLD,"read end %s\n",readitem);
          //printf("OK!!!\n");
          //free(LinVarsele);
        }
        free(sum_cof);
        free(ha_sumele);
      }
      //if(ha_eqint[eqindx]) {
      //  for(lj=0; lj<arSet[ha_eqtime[eqindx]].SetSize; lj++)counteq1[ha_setele[arSet[ha_eqtime[eqindx]].SetBegAdd+lj].setsh]+=nloops/arSet[ha_eqtime[eqindx]].SetSize;
      //} else
      matrow+=nloops;
      eqindx++;

      free(ha_calvar);
      free(LinVars);
      free(arSet);
    }
  }
  //free(counteq1);
  fclose(filehandle);
  return 1;
}

int hlin_dsum(char *formulain, char *commsyntax, hcge_sumcof *sum_cof,ha_cgesetindx *arSet,ha_cgeset *ha_set,uvdim nset,uvdim fdim,int j) {
  char *readitem,*p,*p1,*p2,interchar2[NAMESIZE],argu[TABREADLINE],tempname[NAMESIZE];//,line5[TABREADLINE]
  char interchar[NAMESIZE],interchar1[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],line3[TABREADLINE],line4[TABREADLINE];
  uvdim l,l1,l2,l3,l4,l5,l6,l7;
  int ncur=0,ncuri,i=0,k=0,k1=0,length;
  length=strlen(formulain);
  readitem=formulain;
  //printf("readitem %s\n",readitem);
  //printf("readitem %s\n",commsyntax);
  while (i<length) {
    k=ha_cgefind(readitem,commsyntax);
    //printf("k %d l %d\n",k,length);
    //printf("readitem %s\n",readitem);
    if (k==-1) {
      return 0;
    }
    if (k==0) {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      strcpy(line1,readitem);
      ha_cgecutsum(line);
      //printf("%s\n",line);
      k1=ha_cgefind(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      }
      else {
        if(ha_cgefind(line,",p_")>-1||ha_cgefind(line,"*p_")>-1||ha_cgefind(line,"+p_")>-1||ha_cgefind(line,"-p_")>-1||ha_cgefind(line,"(p_")>-1) {
          i=i+k+4;
          readitem=formulain+i;
        }
        else {
          strcpy(line2,line);
          strcpy(line4,line);
          sprintf(interchar1, "%d", j);
          interchar[0]='\0';
          if (j<10) {
            strcat(interchar,"ha_cgesum000");
          }
          if (9<j&&j<100) {
            strcat(interchar,"ha_cgesum00");
          }
          if (99<j&&j<1000) {
            strcat(interchar,"ha_cgesum0");
          }
          if (999<j&&j<10000) {
            strcat(interchar,"ha_cgesum");
          }
          if (i>10000) {
            strcat(interchar,"ha_cgesum");
            printf("Error: Too many sum\n");
          }
          strcat(interchar,interchar1);
          //ha_cgefrstr(formulain,line,interchar);//good
          //printf("readitem %s\n",formulain);
          strcpy(sum_cof[j].sumname,interchar);
          strcat(interchar,"{");
          strcpy(interchar2,interchar);
          strcpy(line,line1);
          p = strtok(line,"(");
          p = strtok(NULL,",");
          strcpy(sum_cof[j].sumindx,p);
          p = strtok(NULL,",");
          //strcpy(sum_cof[j].sumset,p);
          //printf("sunset %s\n",p);
          for (l7=0; l7<nset; l7++) if(strcmp(p,ha_set[l7].setname)==0) {
              sum_cof[j].sumsetid=l7;
              //sum_cof[j].sumsize=ha_set[l7].size;
              //sum_cof[j].ssetbegadd=ha_set[l7].begadd;
              //sum_cof[j].sumsubsetid=ha_set[l7].subsetid;
              //sum_cof[j].sumsupsetsize=ha_set[l7].supersetsize;
              break;
            }
          //printf("sunsize %d\n",sum_cof[j].sumsize);

          ncur=ha_cgenfind(line2, "{");
          //p = &line2[0];
          l2=0;
          l3=0;
          for (ncuri=0; ncuri<ncur; ncuri++) {
            //strcpy(interchar,interchar2);
            strcpy(line2,line4);
            p = &line2[l2];
            p = strtok(p,"{");
            p = strtok(NULL,"}");
            l2=p-line2;
            strcpy(argu,p);
            strcat(argu,",");
            l=ha_cgenfind(argu, ",");
            //printf("indx %s p %s l %d\n",sum_cof[j].sumindx,p,l);
            if (l<2) {
              if(strcmp(p,sum_cof[j].sumindx)!=0) {
                strcat(interchar,sum_cof[j].sumindx);
                for (l4=0; l4<l3; l4++) {
                  if(strcmp(p,sum_cof[j].dimnames[l4])==0) {
                    break;
                  }
                }
                if (l4==l3||l3==0) {
                  strcpy(sum_cof[j].dimnames[l3],p);
                  l6=0;
                  for (l5=0; l5<fdim-1; l5++) if(strcmp(p,arSet[l5].arIndx)==0) {
                      sum_cof[j].setid[l3]=l5;
                      //strcpy(sum_cof[j].dimsets[l3],arSet[l5].arSet);
                      //sum_cof[j].dims[l3]=arSet[l5].SetSize;
                      //sum_cof[j].subsetid[l3]=arSet[l5].subsetid;
                      //sum_cof[j].supsetsize[l3]=arSet[l5].SuperSetSize;
                      //sum_cof[j].dimssetbegadd[l3]=arSet[l5].SetBegAdd;
                      l6++;
                    }
                  if (l6==0) {
                    interchar1[0]='\0';
                    strcat(interchar1,"sum(");
                    strcat(interchar1,p);
                    strcpy(line3,formulain);
                    line3[readitem-formulain]='\0';
                    l7=ha_cgerevfind(line3,interchar1);
                    p1=&line3[l7+2];
                    p1 = strtok(p1,",");
                    //strcpy(sum_cof[j].dimsets[l3],p1);
                    for (l7=0; l7<nset; l7++) if(strcmp(p1,ha_set[l7].setname)==0) {
                        sum_cof[j].setid[l3]=l7;
                        //sum_cof[j].dims[l3]=ha_set[l7].size;
                        //sum_cof[j].subsetid[l3]=ha_set[l7].subsetid;
                        //sum_cof[j].supsetsize[l3]=ha_set[l7].supersetsize;
                        //sum_cof[j].dimssetbegadd[l3]=ha_set[l7].begadd;
                        break;
                      }
                  }
                  l3++;
                }
              }
            }
            else {
              for (l1=0; l1<l; l1++) {
                if(l1==0) {
                  p = strtok(argu,",");
                }
                else {
                  p = strtok(NULL,",");
                }
                if(strcmp(p,sum_cof[j].sumindx)!=0) {
                  for (l4=0; l4<l3; l4++) {
                    if(strcmp(p,sum_cof[j].dimnames[l4])==0) {
                      break;
                    }
                  }
                  if (l4==l3||l3==0) {
                    strcpy(sum_cof[j].dimnames[l3],p);
                    //printf("p %s\n",p);
                    strcat(interchar,sum_cof[j].dimnames[l3]);
                    l6=0;
                    //printf("fdim %d\n",fdim);
                    for (l5=0; l5<fdim-1; l5++) if(strcmp(p,arSet[l5].arIndx)==0) {
                        sum_cof[j].setid[l3]=l5;
                        //strcpy(sum_cof[j].dimsets[l3],arSet[l5].arSet);
                        //printf("set %s\n",arSet[l5].arSet);
                        //sum_cof[j].dims[l3]=arSet[l5].SetSize;
                        //sum_cof[j].subsetid[l3]=arSet[l5].subsetid;
                        //sum_cof[j].supsetsize[l3]=arSet[l5].SuperSetSize;
                        //sum_cof[j].dimssetbegadd[l3]=arSet[l5].SetBegAdd;
                        l6++;
                        break;
                      }
                    if (l6==0) {
                      interchar1[0]='\0';
                      strcat(interchar1,"sum(");
                      strcat(interchar1,p);
                      strcpy(line3,formulain);
                      line3[readitem-formulain]='\0';
                      l7=ha_cgerevfind(line3,interchar1);
                      p1=&line3[l7+2];
                      p2=strchr(p1,',');
                      //p1 = strtok(p1,",");
                      strncpy(tempname,p1,p2-p1);
                      tempname[p2-p1]='\0';
                      //printf("setname %s\n",sum_cof[j].dimsets[l3]);
                      for (l7=0; l7<nset; l7++) if(strcmp(tempname,ha_set[l7].setname)==0) {
                          sum_cof[j].setid[l3]=l7;
                          //sum_cof[j].dims[l3]=ha_set[l7].size;
                          //sum_cof[j].subsetid[l3]=ha_set[l7].subsetid;
                          //sum_cof[j].supsetsize[l3]=ha_set[l7].supersetsize;
                          //sum_cof[j].dimssetbegadd[l3]=ha_set[l7].begadd;
                          break;
                        }
                    }
                    l3++;
                  }
                }
              }
            }
          }
          if (interchar[strlen(interchar)-1]==',') {
            interchar[strlen(interchar)-1]='}';
          }
          else {
            if (interchar[strlen(interchar)-1]=='{') {
              interchar[strlen(interchar)-1]='\0';
            }
            else {
              strcat(interchar,"}");
            }
          }
          while(ha_cgefrstr(formulain,line4,interchar));
          //printf("readitem %s\n",formulain);
          //printf("line4 %s inter %s\n",line4,interchar);
          sum_cof[j].size=l3;
          return 1;
        }
      }

    }
    else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      strcpy(line1,readitem);
      ha_cgecutsum(line);
      //printf("readitem %s\n",line);
      k1=ha_cgefind(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      }
      else {
        if(ha_cgefind(line,",p_")>-1||ha_cgefind(line,"*p_")>-1||ha_cgefind(line,"+p_")>-1||ha_cgefind(line,"-p_")>-1||ha_cgefind(line,"(p_")>-1) {
          i=i+k+4;
          readitem=formulain+i;
        }
        else {
          strcpy(line2,line);
          strcpy(line4,line);
          //strcpy(line5,line);
          sprintf(interchar1, "%d", j);
          interchar[0]='\0';
          if (j<10) {
            strcat(interchar,"ha_cgesum000");
          }
          if (9<j&&j<100) {
            strcat(interchar,"ha_cgesum00");
          }
          if (99<j&&j<1000) {
            strcat(interchar,"ha_cgesum0");
          }
          if (999<j&&j<10000) {
            strcat(interchar,"ha_cgesum");
          }
          if (i>10000) {
            strcat(interchar,"ha_cgesum");
            printf("Error: Too many sum\n");
          }
          strcat(interchar,interchar1);
          //ha_cgefrstr(formulain,line,interchar);
          //printf("readitem %s\n",readitem);
          strcpy(sum_cof[j].sumname,interchar);
          strcat(interchar,"{");
          strcpy(interchar2,interchar);
          strcpy(line,line1);
          p = strtok(line,"(");
          p = strtok(NULL,",");
          strcpy(sum_cof[j].sumindx,p);
          p = strtok(NULL,",");
          //strcpy(sum_cof[j].sumset,p);
          //printf("sunset %s\n",p);
          for (l7=0; l7<nset; l7++) if(strcmp(p,ha_set[l7].setname)==0) {
              sum_cof[j].sumsetid=l7;
              //sum_cof[j].sumsize=ha_set[l7].size;
              //sum_cof[j].ssetbegadd=ha_set[l7].begadd;
              //sum_cof[j].sumsubsetid=ha_set[l7].subsetid;
              //sum_cof[j].sumsupsetsize=ha_set[l7].supersetsize;
              break;
            }
          //printf("sunsize %d\n",sum_cof[j].sumsize);

          ncur=ha_cgenfind(line2, "{");
          //printf("line2 %s\n",line2);
          //p = &line2[0];
          l2=0;
          l3=0;
          for (ncuri=0; ncuri<ncur; ncuri++) {
            //printf("ncuri %d\n",ncuri);
            strcpy(line2,line4);
            //strcpy(interchar,interchar2);
            p = &line2[l2];
            //printf("indx %s p %s l %d ncuri %d\n",sum_cof[j].sumindx,p,l,ncuri);
            p = strtok(p,"{");
            //printf("line2 %s\n line4 %s\n",line2,line4);
            p = strtok(NULL,"}");
            l2=p-line2;
            strcpy(argu,p);
            strcat(argu,",");
            l=ha_cgenfind(argu, ",");
            //printf("ncuri %d ncur %d argu %s\n",ncuri,ncur,argu);
            if (l<2) {
              if(strcmp(p,sum_cof[j].sumindx)!=0) {
                for (l4=0; l4<l3; l4++) {
                  if(strcmp(p,sum_cof[j].dimnames[l4])==0) {
                    break;
                  }
                }
                if (l4==l3||l3==0) {
                  strcpy(sum_cof[j].dimnames[l3],p);
                  strcat(interchar,sum_cof[j].dimnames[l3]);
                  l6=0;
                  for (l5=0; l5<fdim-1; l5++) if(strcmp(p,arSet[l5].arIndx)==0) {
                      sum_cof[j].setid[l3]=l5;
                      //strcpy(sum_cof[j].dimsets[l3],arSet[l5].arSet);
                      //sum_cof[j].dims[l3]=arSet[l5].SetSize;
                      //sum_cof[j].subsetid[l3]=arSet[l5].subsetid;
                      //sum_cof[j].supsetsize[l3]=arSet[l5].SuperSetSize;
                      //sum_cof[j].dimssetbegadd[l3]=arSet[l5].SetBegAdd;
                      l6++;
                    }
                  if (l6==0) {
                    interchar1[0]='\0';
                    strcat(interchar1,"sum(");
                    strcat(interchar1,p);
                    strcpy(line3,formulain);
                    line3[readitem-formulain]='\0';
                    l7=ha_cgerevfind(line3,interchar1);
                    p1=&line3[l7+2];
                    p1 = strtok(p1,",");
                    //strcpy(sum_cof[j].dimsets[l3],p1);
                    for (l7=0; l7<nset; l7++) if(strcmp(p1,ha_set[l7].setname)==0) {
                        sum_cof[j].setid[l3]=l7;
                        //sum_cof[j].dims[l3]=ha_set[l7].size;
                        //sum_cof[j].subsetid[l3]=ha_set[l7].subsetid;
                        //sum_cof[j].supsetsize[l3]=ha_set[l7].supersetsize;
                        //sum_cof[j].dimssetbegadd[l3]=ha_set[l7].begadd;
                        break;
                      }
                  }
                  l3++;
                }
              }
            }
            else {
              for (l1=0; l1<l; l1++) {
                if(l1==0) {
                  p = strtok(argu,",");
                }
                else {
                  p = strtok(NULL,",");
                }
                //printf("p %s\n",p);
                if(strcmp(p,sum_cof[j].sumindx)!=0) {
                  for (l4=0; l4<l3; l4++) {
                    if(strcmp(p,sum_cof[j].dimnames[l4])==0) {
                      break;
                    }
                  }
                  if (l4==l3||l3==0) {
                    strcpy(sum_cof[j].dimnames[l3],p);
                    //printf("p %s sumindx %s\n",p,sum_cof[j].sumindx);
                    strcat(interchar,sum_cof[j].dimnames[l3]);
                    strcat(interchar,",");
                    l6=0;
                    for (l5=0; l5<fdim-1; l5++) if(strcmp(p,arSet[l5].arIndx)==0) {
                        sum_cof[j].setid[l3]=arSet[l5].setid;//l5;
                        //printf("sum %s dim %s j %d\n",sum_cof[j].sumname,ha_set[sum_cof[j].setid[l3]].setname,j);
                        //strcpy(sum_cof[j].dimsets[l3],arSet[l5].arSet);
                        //sum_cof[j].dims[l3]=arSet[l5].SetSize;
                        //sum_cof[j].subsetid[l3]=arSet[l5].subsetid;
                        //sum_cof[j].supsetsize[l3]=arSet[l5].SuperSetSize;
                        //sum_cof[j].dimssetbegadd[l3]=arSet[l5].SetBegAdd;
                        l6++;
                      }
                    if (l6==0) {
                      interchar1[0]='\0';
                      strcat(interchar1,"sum(");
                      strcat(interchar1,p);
                      strcpy(line3,formulain);
                      line3[readitem-formulain]='\0';
                      l7=ha_cgerevfind(line3,interchar1);
                      p1=&line3[l7+2];
                      p2=strchr(p1,',');
                      //p1 = strtok(p1,",");
                      strncpy(tempname,p1,p2-p1);
                      tempname[p2-p1]='\0';
                      for (l7=0; l7<nset; l7++) if(strcmp(tempname,ha_set[l7].setname)==0) {
                          sum_cof[j].setid[l3]=l7;
                          //sum_cof[j].dims[l3]=ha_set[l7].size;
                          //sum_cof[j].subsetid[l3]=ha_set[l7].subsetid;
                          //sum_cof[j].supsetsize[l3]=ha_set[l7].supersetsize;
                          //sum_cof[j].dimssetbegadd[l3]=ha_set[l7].begadd;
                          break;
                        }
                      //printf("p1 %s dimset %s\n",p1,sum_cof[j].dimsets[l3]);
                    }
                    l3++;
                  }
                }
              }
            }
          }
          if (interchar[strlen(interchar)-1]==',') {
            interchar[strlen(interchar)-1]='}';
          }
          else {
            if (interchar[strlen(interchar)-1]=='{') {
              interchar[strlen(interchar)-1]='\0';
            }
            else {
              strcat(interchar,"}");
            }
          }
          while(ha_cgefrstr(formulain,line4,interchar));
          //printf("readitem %s\n",formulain);
          //printf("line4 %s inter %s\n",line4,interchar);
          sum_cof[j].size=l3;
          return 1;
        }
      }

    }
    else {
      i=i+k+4;
      readitem=formulain+i;
    }
  }
  return 0;
}

int hlin_calsum(char *formulain, char *commsyntax,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele,ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar,hcge_sumcof *sum_cof,int totalsum,ha_cgesumele *ha_sumele,uvadd nsumele,hcge_calvars *ha_calvar,ha_cgesetindx *arSet1,uvdim fdim,int *sumindx,int j, ha_cgetype zerodivide) {
  char *readitem,*p;//,*p1,interchar2[NAMESIZE],line5[TABREADLINE];
  char interchar[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE];//,line3[TABREADLINE],line4[TABREADLINE];//,interchar1[NAMESIZE]
  int i=0,k=0,k1=0,length;//,simpl=0;//,ncur=0,ncuri,l3,l4,l5,l6,l7
  uvdim dcount,setsh,fdimsumcof,l;
  uvadd dcountdim1[4*MAXVARDIM],nloops,l1,l2,l3;
  int ha_calvarsize;
  ha_cgetype vval;
  ha_cgesetindx *arSet2=NULL;
  hcge_calvars *ha_calvar1= NULL;
  uvadd arsetsize;
  length=strlen(formulain);
  readitem=formulain;
  while (i<length) {
    k=ha_cgefind(readitem,commsyntax);
    if (k==-1) {
      return 0;
    }
    if (k==0) {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      ha_cgecutsum(line);
      //printf("here %s\n",line);
      k1=ha_cgefind(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      }
      else {
        if(ha_cgefind(line,",p_")>-1||ha_cgefind(line,"*p_")>-1||ha_cgefind(line,"+p_")>-1||ha_cgefind(line,"-p_")>-1||ha_cgefind(line,"(p_")>-1) {
          i=i+k+4;
          readitem=formulain+i;
        }
        else {
          strcpy(line1,line);
          p=strtok(line,",");
          p=strtok(NULL,",");
          p=strtok(NULL,"\0");
          p[strlen(p)-1]='\0';
          strcpy(line2,p);
          arsetsize=sum_cof[j].size+1;
          ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (arsetsize,sizeof(ha_cgesetindx));
          for (l=0; l<sum_cof[j].size; l++) {
            arSet[l].setid=sum_cof[j].setid[l];
            //strcpy(arSet[l].arSet,sum_cof[j].dimsets[l]);
            strcpy(arSet[l].arIndx,sum_cof[j].dimnames[l]);
            //arSet[l].SetSize=sum_cof[j].dims[l];
            //arSet[l].subsetid=sum_cof[j].subsetid[l];
            //arSet[l].SuperSetSize=sum_cof[j].supsetsize[l];
            //arSet[l].SetBegAdd=sum_cof[j].dimssetbegadd[l];
          }
          nloops=1;
          for (l=0; l<sum_cof[j].size; l++) {
            nloops=nloops*ha_set[arSet[l].setid].size;
            dcount=sum_cof[j].size-l;
            if(dcount==sum_cof[j].size) {
              dcountdim1[dcount-1]=1;
            }
            else {
              dcountdim1[dcount-1]=dcountdim1[dcount]*ha_set[arSet[dcount].setid].size;
            }
          }
          arSet[sum_cof[j].size].setid=sum_cof[j].sumsetid;
          //strcpy(arSet[sum_cof[j].size].arSet,sum_cof[j].sumset);
          strcpy(arSet[sum_cof[j].size].arIndx,sum_cof[j].sumindx);
          //arSet[sum_cof[j].size].SetSize=sum_cof[j].sumsize;
          //arSet[sum_cof[j].size].subsetid=sum_cof[j].sumsubsetid;
          //arSet[sum_cof[j].size].SuperSetSize=sum_cof[j].sumsupsetsize;
          //arSet[sum_cof[j].size].SetBegAdd=sum_cof[j].ssetbegadd;
          fdimsumcof=sum_cof[j].size+1;
          ha_calvarsize=0;
          //printf("p %s\n",p);
          ha_newfparse(p,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,&ha_calvarsize,arSet,fdimsumcof);
          //for(l=0; l<ha_calvarsize; l++) for(l2=0; l2<10; l2++) {
          //printf("name %s oper %d tmptype %d begadd %d ad %d\n",ha_calvar[l].TmpVarName,ha_calvar[l].Oper,ha_calvar[l].Var1Type,ha_calvar[l].Var1BegAdd,ha_calvar[l].Var1ADims[l2]);
          //}
          //for(l=0;l<ha_calvarsize;l++) for(l2=0;l2<10;l2++) printf("tmptype %d begadd %d ad %d\n",ha_calvar[l].TmpVarType,ha_calvar[l].Var1BegAdd,ha_calvar[l].Var1ADims[l2]);
        #pragma omp parallel private(l3,l1,l2,dcount,setsh,vval,arSet2,ha_calvar1) shared(ha_cofvar,arSet,ha_sumele)
        {
        if(omp_get_thread_num()!=0){
          arSet2=realloc(arSet2,arsetsize*sizeof(ha_cgesetindx));
          memcpy(arSet2,arSet,arsetsize*sizeof(ha_cgesetindx));
          ha_calvar1=realloc(ha_calvar1,ha_calvarsize*sizeof(hcge_calvars));
          memcpy(ha_calvar1,ha_calvar,ha_calvarsize*sizeof(hcge_calvars));
        }else{
          ha_calvar1=ha_calvar;
          arSet2=arSet;
        }
        #pragma omp for
          for (l3=0; l3<nloops; l3++) {
            l2=l3;
            for (dcount=0; dcount<sum_cof[j].size; dcount++) {
              setsh=(uvadd) l2/dcountdim1[dcount];
              arSet2[dcount].indx=setsh;
              l2=l2-setsh*dcountdim1[dcount];
            }
            vval=0;
            for (l1=0; l1<ha_set[sum_cof[j].sumsetid].size; l1++) {
              arSet2[sum_cof[j].size].indx=l1;
              vval+=ha_newfpcal(ha_cofvar,ha_set,ha_setele,ha_sumele,ha_calvar1,ha_calvarsize,arSet2,fdimsumcof,zerodivide);
            }
            ha_sumele[*sumindx+l3].varval=vval;
            //printf("*sumindx+l %d sum %f\n",*sumindx+l,vval);
          }
        if(omp_get_thread_num()!=0){
          free(arSet2);
          arSet2=NULL;
          free(ha_calvar1);
          ha_calvar1=NULL;
        }else{
          ha_calvar1=NULL;
          arSet2=NULL;
        }
        }
          *sumindx=*sumindx+nloops;
          strcpy(interchar,sum_cof[j].sumname);
          strcat(interchar,"{");
          for (l=0; l<sum_cof[j].size; l++) {
            strcat(interchar,sum_cof[j].dimnames[l]);
            strcat(interchar,",");
          }
          if (interchar[strlen(interchar)-1]==',') {
            interchar[strlen(interchar)-1]='}';
          }
          else {
            if (interchar[strlen(interchar)-1]=='{') {
              interchar[strlen(interchar)-1]='\0';
            }
            else {
              strcat(interchar,"}");
            }
          }
          while(ha_cgefrstr(formulain,line1,interchar)!=NULL);
          free(arSet);
          return 1;
        }
      }

    }
    else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      ha_cgecutsum(line);
      //printf("sum %s\n",line);//de check
      k1=ha_cgefind(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      }
      else {
        if(ha_cgefind(line,",p_")>-1||ha_cgefind(line,"*p_")>-1||ha_cgefind(line,"+p_")>-1||ha_cgefind(line,"-p_")>-1||ha_cgefind(line,"(p_")>-1) {
          i=i+k+4;
          //printf("i %d\n",i);
          readitem=formulain+i;
        }
        else {
          strcpy(line1,line);
          p=strtok(line,",");
          p=strtok(NULL,",");
          p=strtok(NULL,"\0");
          p[strlen(p)-1]='\0';
          //strcpy(line2,p);
          //printf("p %s\n",p);
          //if (strpbrk(p,"*+-^/)")==NULL) {
          //simpl=1;
          //}
          arsetsize=sum_cof[j].size+1;
          ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (arsetsize,sizeof(ha_cgesetindx));
          for (l=0; l<sum_cof[j].size; l++) {
            arSet[l].setid=sum_cof[j].setid[l];
            //strcpy(arSet[l].arSet,sum_cof[j].dimsets[l]);
            strcpy(arSet[l].arIndx,sum_cof[j].dimnames[l]);
            //arSet[l].subsetid=sum_cof[j].subsetid[l];
            //arSet[l].SuperSetSize=sum_cof[j].supsetsize[l];
            //arSet[l].SetBegAdd=sum_cof[j].dimssetbegadd[l];
            //arSet[l].SetSize=sum_cof[j].dims[l];
          }
          nloops=1;
          for (l=0; l<sum_cof[j].size; l++) {
            nloops=nloops*ha_set[arSet[l].setid].size;//sum_cof[j].dims[l];
            dcount=sum_cof[j].size-l;
            if(dcount==sum_cof[j].size) {
              dcountdim1[dcount-1]=1;
            }
            else {
              dcountdim1[dcount-1]=dcountdim1[dcount]*ha_set[arSet[dcount].setid].size;
            }
          }
          arSet[sum_cof[j].size].setid=sum_cof[j].sumsetid;
          //strcpy(arSet[sum_cof[j].size].arSet,sum_cof[j].sumset);
          strcpy(arSet[sum_cof[j].size].arIndx,sum_cof[j].sumindx);
          //arSet[sum_cof[j].size].SetSize=sum_cof[j].sumsize;
          //arSet[sum_cof[j].size].subsetid=sum_cof[j].sumsubsetid;
          //arSet[sum_cof[j].size].SuperSetSize=sum_cof[j].sumsupsetsize;
          //arSet[sum_cof[j].size].SetBegAdd=sum_cof[j].ssetbegadd;
          fdimsumcof=sum_cof[j].size+1;
          ha_calvarsize=0;
          //printf("p %s nloop %d\n",p,nloops);
          ha_newfparse(p,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,&ha_calvarsize,arSet,fdimsumcof);
          //printf("p %s nloop %d\n",p,nloops);
          //for(l=0; l<ha_calvarsize; l++) for(l2=0; l2<10; l2++) {
          //printf("name %s oper %d tmptype %d begadd %d ad %d\n",ha_calvar[l].TmpVarName,ha_calvar[l].Oper,ha_calvar[l].Var1Type,ha_calvar[l].Var1BegAdd,ha_calvar[l].Var1ADims[l2]);
          //}
          //for(l=0;l<ha_calvarsize;l++) for(l2=0;l2<10;l2++) printf("tmptype %d begadd %d ad %d\n",ha_calvar[l].TmpVarType,ha_calvar[l].Var1BegAdd,ha_calvar[l].Var1ADims[l2]);
        #pragma omp parallel private(l3,l1,l2,dcount,setsh,vval,arSet2,ha_calvar1) shared(ha_cofvar,arSet,ha_sumele)
        {
        if(omp_get_thread_num()!=0){
          arSet2=realloc(arSet2,arsetsize*sizeof(ha_cgesetindx));
          memcpy(arSet2,arSet,arsetsize*sizeof(ha_cgesetindx));
          ha_calvar1=realloc(ha_calvar1,ha_calvarsize*sizeof(hcge_calvars));
          memcpy(ha_calvar1,ha_calvar,ha_calvarsize*sizeof(hcge_calvars));
        }else{
          ha_calvar1=ha_calvar;
          arSet2=arSet;
        }
        #pragma omp for
          for (l3=0; l3<nloops; l3++) {
            l2=l3;
            for (dcount=0; dcount<sum_cof[j].size; dcount++) {
              setsh=(uvadd) l2/dcountdim1[dcount];
              arSet2[dcount].indx=setsh;
              l2=l2-setsh*dcountdim1[dcount];
            }
            vval=0;
            //printf("l3 %d l2 %d sh %d\n",l3,l2,setsh);
            for (l1=0; l1<ha_set[sum_cof[j].sumsetid].size; l1++) {
              //printf("l1 %d\n",l1);
              arSet2[sum_cof[j].size].indx=l1;
              //strcpy(p,line2);
              vval+=ha_newfpcal(ha_cofvar,ha_set,ha_setele,ha_sumele,ha_calvar1,ha_calvarsize,arSet2,fdimsumcof,zerodivide);
              //if (simpl==1) {
              //vval=vval+hcge_varrepl(p,ha_cofvar,ncofele,ha_cof,ncof,ha_var,nvar,ha_setele,sum_cof,totalsum,ha_sumele,nsumele,arSet,fdimsumcof);//hcge_ssimplrpl(p,ha_cofvar,ncofele,ha_cof,ncof,ha_var,nvar,sum_cof,totalsum,ha_sumele,nsumele,arSet,fdimsumcof);
              //} else {
              //vval=vval+ha_cgefparse(p,ha_cofvar,ncofvar,ncofele,ha_cof,ncof,ha_var,nvar,ha_setele,sum_cof,totalsum,ha_sumele,nsumele,vartemppow,vartempmuldiv,vartempplu,arSet,fdimsumcof,zerodivide);
              //}
            }
            //printf("*sumindx+l %d sum %f sumele %d hasum %f\n",*sumindx+l3,vval,nsumele,ha_sumele[0].varval);
            ha_sumele[(uvadd)*sumindx+l3].varval=vval;//ha_sumele[*sumindx+l2].varval=vval;
          }
        if(omp_get_thread_num()!=0){
          free(arSet2);
          arSet2=NULL;
          free(ha_calvar1);
          ha_calvar1=NULL;
        }else{
          ha_calvar1=NULL;
          arSet2=NULL;
        }
        }
          //printf("OK!!!\n");
          *sumindx=*sumindx+nloops;
          strcpy(interchar,sum_cof[j].sumname);
          strcat(interchar,"{");
          for (l=0; l<sum_cof[j].size; l++) {
            strcat(interchar,sum_cof[j].dimnames[l]);
            strcat(interchar,",");
          }
          if (interchar[strlen(interchar)-1]==',') {
            interchar[strlen(interchar)-1]='}';
          }
          else {
            if (interchar[strlen(interchar)-1]=='{') {
              interchar[strlen(interchar)-1]='\0';
            }
            else {
              strcat(interchar,"}");
            }
          }
          while(ha_cgefrstr(formulain,line1,interchar));
          free(arSet);
          return 1;
        }
      }
    }
    else {
      i=i+k+4;
      readitem=formulain+i;
    }
  }
  return 0;
}

int hlin_replsum(char *formulain, char *commsyntax,int LinIndx, hcge_linvars *LinVars,hcge_cof *ha_var) {
  char *readitem,line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],line3[TABREADLINE],argu[TABREADLINE],*p0,*p1;
  int i=0,k=0,k1=0,k2,k3,length,tsum,j,l;
  bool IsRemSum;
  length=strlen(formulain);
  readitem=formulain;
  //printf("readitem %s\n",readitem);
  //printf("readitem %s\n",commsyntax);
  while (i<length) {
    k=ha_cgefind(formulain+i,commsyntax);
    //printf("k %d l %d\n",k,length);
    //printf("readitem %s\n",readitem);
    if (k==-1) {
      return 0;
    }
    if (k==0) {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      ha_cgecutsum(line);
      strcpy(line3,line);
      if(ha_cgefind(line,",p_")>-1||ha_cgefind(line,"*p_")>-1||ha_cgefind(line,"+p_")>-1||ha_cgefind(line,"-p_")>-1||ha_cgefind(line,"(p_")>-1) {
        tsum=hcge_nsum(line,commsyntax);
        IsRemSum=false;
        for (j=tsum; j>1; j--) {
          k2=1;
          k1=0;
          while(k2<j) {
            k3=k1;
            k1=ha_cgefind(line+k1+3,commsyntax);
            k1=k3+k1+3;//printf("k1 %d line %s\n",k1,line+k1-1);
            //printf("k1 %s\n",line+k1);
            if (line[k1-1]=='+'||line[k1-1]=='-'||line[k1-1]=='*'||line[k1-1]=='/'||line[k1-1]=='^'||line[k1-1]=='('||line[k1-1]==',') k2++;//printf("k1 %d line %s\n",k1,line+k1-1);}
          }
          strcpy(line1,&line[k1]);
          ha_cgecutsum(line1);
          p0=strchr(line1,'(');
          p0++;
          p1=strchr(line1,',');
          strncpy(argu,p0,p1-p0);
          argu[p1-p0]='\0';
          for(l=0; l<ha_var[LinVars[LinIndx].LinVarIndx].size; l++) {
            if(strcmp(argu,LinVars[LinIndx].dimnames[l])==0) {
              p1++;
              p1=strchr(p1,',');
              p1++;
              strncpy(line2,p1,strlen(line1)-(p1-line1)-1);
              line2[strlen(line1)-(p1-line1)-1]='\0';
              ha_cgefrstr(line,line1,line2);
              IsRemSum=true;
              break;
            }
          }
        }
        p0=strchr(line,'(');
        p0++;
        p1=strchr(line,',');
        strncpy(argu,p0,p1-p0);
        argu[p1-p0]='\0';
        //printf("argu %s\n",argu);
        for(l=0; l<ha_var[LinVars[LinIndx].LinVarIndx].size; l++) {
          //printf("argu %s var %s\n",argu,LinVars[LinIndx].dimnames[l]);
          if(strcmp(argu,LinVars[LinIndx].dimnames[l])==0) {
            p1++;
            p1=strchr(p1,',');
            p1++;
            line2[0]='(';
            strncpy(&line2[1],p1,strlen(line)-(p1-line)-1);
            line2[strlen(line)-(p1-line)]=')';
            line2[strlen(line)-(p1-line)+1]='\0';
            strcpy(line,line2);
            IsRemSum=true;
            break;
          }
        }
        if(IsRemSum) {
          ha_cgefrstr(formulain,line3,line);
          i=i+k+strlen(line);
        }
        else i=i+k+strlen(line3);
        //return 1;
      }
      else {
        ha_cgefrstr(formulain,line,"0");
        i=i+k+1;
        //return 1;
      }
    }
    else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      ha_cgecutsum(line);
      strcpy(line3,line);
      if(ha_cgefind(line,",p_")>-1||ha_cgefind(line,"*p_")>-1||ha_cgefind(line,"+p_")>-1||ha_cgefind(line,"-p_")>-1||ha_cgefind(line,"(p_")>-1) {
        tsum=hcge_nsum(line,commsyntax);
        IsRemSum=false;
        for (j=tsum; j>1; j--) {
          k2=1;
          k1=0;
          while(k2<j) {
            k3=k1;
            k1=ha_cgefind(line+k1+3,commsyntax);
            k1=k3+k1+3;//printf("k1 %d line %s\n",k1,line+k1-1);
            if (line[k1-1]=='+'||line[k1-1]=='-'||line[k1-1]=='*'||line[k1-1]=='/'||line[k1-1]=='^'||line[k1-1]=='('||line[k1-1]==',') k2++;//printf("k1 %d line %s\n",k1,line+k1-1);}
          }
          strcpy(line1,&line[k1]);
          //printf("line1 %s\n",line1);
          ha_cgecutsum(line1);
          //printf("line1 %s\n",line1);
          p0=strchr(line1,'(');
          p0++;
          p1=strchr(line1,',');
          strncpy(argu,p0,p1-p0);
          argu[p1-p0]='\0';
          //printf("line1 %s\n",line1);
          for(l=0; l<ha_var[LinVars[LinIndx].LinVarIndx].size; l++) {
            if(strcmp(argu,LinVars[LinIndx].dimnames[l])==0) {
              p1++;
              p1=strchr(p1,',');
              p1++;
              strncpy(line2,p1,strlen(line1)-(p1-line1)-1);
              line2[strlen(line1)-(p1-line1)-1]='\0';
              //printf("line %s\n line1 %s\n line2 %s\n",line,line1,line2);
              ha_cgefrstr(line,line1,line2);
              IsRemSum=true;
              break;
            }
          }
        }
        //printf("line %s\n",line);
        p0=strchr(line,'(');
        p0++;
        p1=strchr(line,',');
        strncpy(argu,p0,p1-p0);
        argu[p1-p0]='\0';
        for(l=0; l<ha_var[LinVars[LinIndx].LinVarIndx].size; l++) {
          if(strcmp(argu,LinVars[LinIndx].dimnames[l])==0) {
            p1++;
            p1=strchr(p1,',');
            p1++;
            line2[0]='(';
            strncpy(&line2[1],p1,strlen(line)-(p1-line)-1);
            line2[strlen(line)-(p1-line)]=')';
            line2[strlen(line)-(p1-line)+1]='\0';
            strcpy(line,line2);
            IsRemSum=true;
            break;
          }
        }
        if(IsRemSum) {
          ha_cgefrstr(formulain,line3,line);
          i=i+k+strlen(line);
        }
        else i=i+k+strlen(line3);
        //return 1;
      }
      else {
        ha_cgefrstr(formulain,line,"0");
        i=i+k+1;
        //return 1;
      }
    }
    else {
      i=i+k+4;
      readitem=formulain+i;
    }
  }
  return 0;
}

int hlin_rlinone(char *formulain,hcge_linvars *LinVars,int linindx,hcge_cof *ha_var) {
  char line[TABREADLINE];
  int i1,i2,i3,l;//,l,d;
  line[0]='p';
  line[1]='_';
  line[2]='\0';
  l=strlen(formulain);
  strcat(line,LinVars[linindx].LinVarName);
  //printf("%d varname %s\n",ha_var[LinVars[linindx].LinVarIndx].size,line);
  i3=0;
  i1=0;
  //printf("OK!!!\n");
  while (i3==0) {
    i1=ha_cgefind(formulain+i1,line);
    if (i1==-1) break;
    if (i1==0) {
      if(ha_var[LinVars[linindx].LinVarIndx].size==0) {
        i2=l-1;
        //printf("%s\n",line);
      }
      else
        i2=ha_cgefind(formulain,"}");
      formulain[0]='1';
      memmove(formulain+1,formulain+i2+1,l-i2);
      i3=1;
    }
    else {
      if (formulain[i1-1]=='+'||formulain[i1-1]=='-'||formulain[i1-1]=='*'||formulain[i1-1]=='/'||formulain[i1-1]=='^'||formulain[i1-1]=='('||formulain[i1-1]==',') {
        if(ha_var[LinVars[linindx].LinVarIndx].size==0) {
          i2=strlen(line)-1;
          //printf("%s\n",line);
        }
        else i2=ha_cgefind(formulain+i1,"}");
        formulain[i1]='1';
        memmove(formulain+i1+1,formulain+i1+i2+1,l-i2-i1);
        i3=1;
      }
      else i1++;
    }
  }
  return 1;
}

int NewMatvarRead(char *fname, char *commsyntax,ha_cgeset *ha_set,uvdim nset,ha_cgesetele *ha_setele,hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar,ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele,ha_cgeexovar *ha_cgeshock,bool *var_inter,hcge_cof *ha_eq,bool *ha_eqint,uvdim *eq_orderintra,uvdim *eq_orderreg,uvadd allregset,uvadd alltimeset,uvdim *orderintra,uvdim *orderreg) {
  FILE * filehandle;
  char tline[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[TABREADLINE],lintmp[TABREADLINE];//,*p1=NULL;
  char *readitem=NULL,*p=NULL,*p1=NULL;//,*p2=NULL,*varpnts;
  ha_cgetype zerodivide=0;
  uvdim fdim,np,i4;
  uvadd rowindx,j,l,l1,lj,dcountdim1[4*MAXVARDIM],dcountdim2[4*MAXVARDIM],nloops,nloopslin,li3,l2,matrow,eqindx=0,matsize,l01,j01;//,sizelinvars,totlinvars,templinvars
  int sj,i,i3,nlinvars,lvar,lvar1,lvar2,lvar3,lvar4,dcount,fdimlin=0,leadlag,varindx1,varindx2,tempint;
  ha_cgetype dimmat[4*MAXVARDIM];
  //bool *varchoose= (bool *) calloc (nvar,sizeof(bool));
  //printf("varchoose %d\n",varchoose[70]);

  filehandle = fopen(fname,"r");
  matrow=0;
  //printf("OK!!!!\n");

  while (ha_cgertabl1(commsyntax,filehandle,line,ha_cofvar,ha_cof,ncof,&zerodivide,TABREADLINE)) {
    //PetscPrintf(PETSC_COMM_WORLD,"line %s\n",line);
    if (strstr(line,"(default")==NULL) {
      ha_cgefrstr1(line, commsyntax, "");
      ha_cgefrstr1(line, "(linear)", "");
      while (ha_cgefrstr(line,"  ", " "));
      while (ha_cgefrchr(line, '[', '('));
      while (ha_cgefrchr(line, ']', ')'));
      while (ha_cgefrchr(line, '{', '('));
      while (ha_cgefrchr(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=ha_cgenfind(line, "(all,");
      //printf("line %s\n",line);
      if (fdim==0) {
        readitem = strtok(line+1," ");
        strcpy(ha_eq[eqindx].cofname,readitem);
        //printf("eq name %s\n",readitem);
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      else {
        //printf("line %s\n",line+1);
        readitem = strtok(line+1,"(");
        //printf("read %s\n",readitem);
        strcpy(ha_eq[eqindx].cofname,readitem);
        //printf("eq name %s\n",readitem);
        strcpy(line,linecopy);
        i=ha_cgerevfind(line, "(all,");
        readitem=line+i;
        readitem = strtok(readitem,")");
        readitem = strtok(NULL,"=");
        //printf("left %s\n",readitem);
        strcpy(vname,readitem);
        //printf("left %s linecopy %s\n",linecopy1,linecopy);
        strcpy(line,linecopy);
        //printf("left %s\n",line);
        readitem = strtok(line,"=");
        //printf("read %s\n",readitem);
        readitem = strtok(NULL,";");
        //printf("left %s\n",readitem);
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      while (ha_cgefrstr(readitem," ", ""));
      while (ha_cgerecovar(readitem)==1);
      hnew_intrpl(readitem);
      //printf("read %s\n",readitem);
      strcpy(tline,readitem);
      strcpy(line1,readitem);
      np=ha_cgenfind(readitem,"p_");
      hcge_linvars *LinVars= (hcge_linvars *) calloc (np+1,sizeof(hcge_linvars));
      i3=0;
      //sizelinvars=0;
      lvar=0;
      for (i=0; i<np; i++) {
        //varindx1=0;
        varindx2=0;
        while(-1<0) {
          varindx1=ha_cgefind(readitem+varindx2,"p_");
          if(varindx1==-1) break;
          varindx2=varindx2+varindx1;
          if(varindx2==0||readitem[varindx2-1]=='*'||readitem[varindx2-1]=='+'||readitem[varindx2-1]=='-'||readitem[varindx2-1]=='('||readitem[varindx2-1]==',') break;
          else varindx2++;
        }
        if(varindx1==-1) break;
        lvar=lvar+varindx2;
        readitem=readitem+varindx2+2;
        //printf("varbeg %s\n",readitem);
        p=strpbrk(readitem,"{+*-/^)");
        if(p!=NULL) {
          strncpy(vname,readitem,p-readitem);
          vname[p-readitem]='\0';
          //printf("vname %s\n",vname);
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(ha_var[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          p=strpbrk(readitem,"}+*-/^)");
          if (*p=='}') {
            p++;
          }
          strncpy(tline,readitem,p-readitem);
          tline[p-readitem]='\0';
          switch (ha_var[l].size) {
          case 0:
            break;
          case 1:
            p = strtok(tline,"{");
            p = strtok(NULL,"}");
            leadlag=0;
            hnew_arset(p,&leadlag);
            if(leadlag!=0)var_inter[l]=true;//printf("var %s\n",ha_var[l].cofname);}
            strcpy(LinVars[i3].dimnames[0],p);
            //LinVars[i3].dimleadlag[0]=leadlag;
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=ha_cgefind(linecopy,lintmp);
            //printf("lin %s\n",lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[0],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              //printf("lin1 %s\n",lintmp);
              lvar1=ha_cgenfind(linecopy,lintmp);
              lvar3=ha_cgefind(linecopy,lintmp);
              //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
                  if (lvar4>-1&&lvar4<lvar) {
                    lvar3=lvar3+lvar4+4;
                  }
                  else {
                    break;
                  }
                }
              //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
              p1=&linecopy[0]+lvar3;
              //p1=&linecopy[0]+ha_cgefind(linecopy,lintmp);
              p1=p1+strlen(lintmp);
              //printf("lin %s\n",lintmp);
              strncpy(LinVars[i3].dimsetnames[0],p1,strchr(p1,',')-p1);
            }
            //printf("vname %s lin %s\n",vname,LinVars[i3].dimsetnames[0]);
            //printf("lin %s indx %s\n",ha_var[l].cofname,p);
            break;
          default:
            p = strtok(tline,"{");
            for (i4=0; i4<ha_var[l].size-1; i4++) {
              p = strtok(NULL,",");
              leadlag=0;
              hnew_arset(p,&leadlag);
              if(leadlag!=0)var_inter[l]=true;//printf("var %s\n",ha_var[l].cofname);}
              strcpy(LinVars[i3].dimnames[i4],p);
              //printf("p %s\n",p);
              //LinVars[i3].dimleadlag[i4]=leadlag;
              strcpy(lintmp,"(all,");
              strcat(lintmp,p);
              strcat(lintmp,",");
              l1=ha_cgefind(linecopy,lintmp);
              //printf("lin %s\n",lintmp);
              if (l1>-1) {
                p1=&linecopy[0]+l1;
                strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
              }
              else {
                strcpy(lintmp,"sum(");
                strcat(lintmp,p);
                strcat(lintmp,",");
                //printf("lin1 %s\n",lintmp);
                lvar1=ha_cgenfind(linecopy,lintmp);
                lvar3=ha_cgefind(linecopy,lintmp);
                //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
                if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                    lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
                    if (lvar4>-1&&lvar4<lvar) {
                      lvar3=lvar3+lvar4+4;
                    }
                    else {
                      break;
                    }
                  }
                //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
                p1=&linecopy[0]+lvar3;
                //p1=&linecopy[0]+ha_cgefind(linecopy,lintmp);
                p1=p1+strlen(lintmp);
                //printf("lin %s p1 %s len %d\n",lintmp,p1,strchr(p1,',')-p1-strlen(lintmp));
                strncpy(LinVars[i3].dimsetnames[i4],p1,strchr(p1,',')-p1);
                //printf("lin %s\n",LinVars[i3].dimsetnames[i4]);
              }
              //printf("vname %s lin %s\n",vname,LinVars[i3].dimsetnames[i4]);
            }
            p = strtok(NULL,"}");
            leadlag=0;
            hnew_arset(p,&leadlag);
            if(leadlag!=0)var_inter[l]=true;//printf("var %s\n",ha_var[l].cofname);}
            strcpy(LinVars[i3].dimnames[i4],p);
            //LinVars[i3].dimleadlag[i4]=leadlag;
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=ha_cgefind(linecopy,lintmp);
            //printf("lin %s\n",lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=ha_cgenfind(linecopy,lintmp);
              lvar3=ha_cgefind(linecopy,lintmp);
              //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
                  if (lvar4>-1&&lvar4<lvar) {
                    lvar3=lvar3+lvar4+4;
                  }
                  else {
                    break;
                  }
                }
              //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
              p1=&linecopy[0]+lvar3;
              p1=p1+strlen(lintmp);
              //printf("lin %s p1 %s\n",lintmp,p1);
              strncpy(LinVars[i3].dimsetnames[i4],p1,strchr(p1,',')-p1);
            }
            //printf("vname %s lin %s l1 %d\n",vname,LinVars[i3].dimsetnames[i4],l1);
            break;
          }
          i3++;
          //}
          //LinVars[i][p-readitem]='\0';
        }
        else {
          strcpy(vname,readitem);
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(ha_var[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          i3++;
          //}
          //break;
        }
      }
      nlinvars=i3;
      for (i=0; i<nlinvars; i++) i3+=ha_var[LinVars[i].LinVarIndx].size;
      ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (fdim+i3+1,sizeof(ha_cgesetindx));
      bool *linvarrcount= (bool *) calloc (nlinvars,sizeof(bool));
      //bool *linvartcount= (bool *) calloc (nlinvars,sizeof(bool));

      strcpy(line,linecopy);
      //printf("line %s\n",line);
      nloops=1;
      if (fdim>0) {
        for (i=0; i<fdim; i++) {
          if(i==0) {
            readitem = strtok(line,",");
          }
          else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
          strcpy(arSet[i].arIndx,readitem);
          readitem = strtok(NULL,")");
          for (i4=0; i4<nset; i4++) if(strcmp(readitem,ha_set[i4].setname)==0) {
              arSet[i].setid=i4;
              ha_eq[eqindx].setid[i]=i4;
              break;
            }
          nloops=nloops*ha_set[arSet[i].setid].size;
        }
      }
      ha_eq[eqindx].size=fdim;
      ha_eq[eqindx].matsize=nloops;
      if(alltimeset>=0) {
        if(eqindx==0) {
          for(i=0; i<nvar; i++) {
            if(var_inter[i]) continue;
            j=0;
            for(i4=0; i4<ha_var[i].size; i4++)if(ha_set[ha_var[i].setid[i4]].intertemp) {
                orderintra[i]=i4;
                j++;
                break;
              }
            if(j==0)var_inter[i]=true;
          }
        }
        j=0;
        for (i=0; i<fdim; i++)if(ha_set[ha_eq[eqindx].setid[i]].intertemp) {
            j++;
            eq_orderintra[eqindx]=i;
          }
        if(j==0) {
          if(allregset<0)for (i4=0; i4<nlinvars; i4++) var_inter[LinVars[i4].LinVarIndx]=true;
          ha_eqint[eqindx]=true;
        }

        /*if(eqindx==0) {
          for(i=0; i<nvar; i++) {
            j=0;
            for(i4=0; i4<ha_var[i].size; i4++)if(ha_set[ha_var[i].setid[i4]].intertemp) {
                orderreg[i]=i4;
                j++;
                break;
              }
            if(j==0)var_inter[i]=true;
          }
        }
        j=0;
        for(i4=0; i4<nlinvars; i4++)for(i=0; i<ha_var[LinVars[i4].LinVarIndx].size; i++) if(ha_set[ha_var[LinVars[i4].LinVarIndx].setid[i]].intertemp)j++;
        if(j==0)ha_eqint[eqindx]=true;
        j=0;
        for (i=0; i<fdim; i++)if(ha_set[ha_eq[eqindx].setid[i]].intertemp) {
            j++;
            eq_orderintra[eqindx]=i;
          }
        if(j==0)ha_eqint[eqindx]=true;
        if(j>=2) {
          eq_orderintra[eqindx]=-1;
          for (i=0; i<fdim; i++) {
            for(i4=0; i4<nlinvars; i4++) {
              for(i3=0; i3<ha_var[LinVars[i4].LinVarIndx].size; i3++) {
                if(strcmp(LinVars[i4].dimnames[i3],arSet[i].arIndx)==0&&ha_set[ha_eq[eqindx].setid[i]].intertemp&&ha_set[ha_var[LinVars[i4].LinVarIndx].setid[i3]].intertemp) {
                  eq_orderintra[eqindx]=i;
                  break;
                }
              }
              if(eq_orderintra[eqindx]>=0)break;
            }
            if(eq_orderintra[eqindx]>=0)break;
          }
          //printf("eq %s eqreg %d\n",ha_eq[eqindx].cofname,eq_orderreg[eqindx]);
        }
        matsize=0;
        for(i4=0; i4<nlinvars; i4++)
          for(i3=0; i3<ha_var[LinVars[i4].LinVarIndx].size; i3++) {
            if(strcmp(LinVars[i4].dimnames[i3],arSet[eq_orderintra[eqindx]].arIndx)!=0&&orderintra[LinVars[i4].LinVarIndx]==i3&&!var_inter[LinVars[i4].LinVarIndx]) {
              matsize+=ha_var[LinVars[i4].LinVarIndx].matsize;
            }
          }
        if(matsize>ha_eq[eqindx].matsize)ha_eqint[eqindx]=true;
        else {
          for(i4=0; i4<nlinvars; i4++)
            for(i3=0; i3<ha_var[LinVars[i4].LinVarIndx].size; i3++) {
              if(strcmp(LinVars[i4].dimnames[i3],arSet[eq_orderintra[eqindx]].arIndx)!=0&&orderintra[LinVars[i4].LinVarIndx]==i3&&!var_inter[LinVars[i4].LinVarIndx]) {
                var_inter[LinVars[i4].LinVarIndx]=true;
              }
            }
        }*/
      }
      //for(i=0;i<39;i++)printf("i %d neq %d\n",i,ha_eqint[i]);
      if(allregset>=0) {
        if(eqindx==0) {
          for(i=0; i<nvar; i++) {
            j=0;
            for(i4=0; i4<ha_var[i].size; i4++)if(ha_set[ha_var[i].setid[i4]].regional) {
                orderreg[i]=i4;
                //if(i==23)printf("!!!!!!!!!!!!!!!!!!!!!!! var %s\n",ha_var[i].cofname);
                //if(i==230)orderreg[i]=2;
                j++;
                //break;
              }
            if(j==0)var_inter[i]=true;
            //if(j==2)printf("!!!!!!!!!!!!!!!!!!!!!!! i %d var %s\n",i,ha_var[i].cofname);
          }
        }
        j=0;
        for(i4=0; i4<nlinvars; i4++)for(i=0; i<ha_var[LinVars[i4].LinVarIndx].size; i++) if(ha_set[ha_var[LinVars[i4].LinVarIndx].setid[i]].regional)j++;
        if(j==0)ha_eqint[eqindx]=true;
        j=0;
        for (i=0; i<fdim; i++)if(ha_set[ha_eq[eqindx].setid[i]].regional) {
            j++;
            eq_orderreg[eqindx]=i;
          }
        //printf("OK!!!!\n");printf("eqindx %d s of %d eqindx %d\n",eqindx,sizeof(ha_eqint),1);
        if(j==0) {
          //for (i4=0; i4<nlinvars; i4++) var_inter[LinVars[i4].LinVarIndx]=true;
          //printf("var %s\n",ha_var[LinVars[i4].LinVarIndx].cofname);
          ha_eqint[eqindx]=true;
        }
        //if(strcmp(ha_eq[eqindx].cofname,"exprices")==0)printf("j %d\n",j);
        if(j>=2) {
          tempint=eq_orderreg[eqindx];
          eq_orderreg[eqindx]=-1;
          for (i=0; i<fdim; i++)dimmat[i]=-1;
          for (i=0; i<fdim; i++) {
            memset(linvarrcount,0,nlinvars*sizeof(bool));
            for(i4=0; i4<nlinvars; i4++) {
              //printf("i4 %d linvar %s\n",i4,LinVars[i4].LinVarName);
              for(lvar4=0; lvar4<i4; lvar4++)if(strcmp(LinVars[i4].LinVarName,LinVars[lvar4].LinVarName)==0)break;//{printf("here!!! linvarrcount[lvar4] %d\n",linvarrcount[lvar4]);
              if(lvar4<i4&&linvarrcount[lvar4]) {
                linvarrcount[i4]=true;
                continue;
              }
              for(i3=0; i3<ha_var[LinVars[i4].LinVarIndx].size; i3++) {
                if(strcmp(LinVars[i4].dimnames[i3],arSet[i].arIndx)==0&&ha_set[ha_eq[eqindx].setid[i]].regional&&orderreg[LinVars[i4].LinVarIndx]==i3) { //ha_set[ha_var[LinVars[i4].LinVarIndx].setid[i3]].regional) {
                  dimmat[i]+=ha_var[LinVars[i4].LinVarIndx].matsize;
                  l01=0;
                  for(j01=0; j01<ha_var[LinVars[i4].LinVarIndx].matsize; j01++)if(ha_cgeshock[ha_var[LinVars[i4].LinVarIndx].begadd+j01].ShockId)l01++;
                  dimmat[i]-=l01;
                  linvarrcount[i4]=true;
                  if(lvar4<i4)linvarrcount[lvar4]=true;
                  //break;
                }
              }
              //if(eq_orderreg[eqindx]>=0)break;
            }
            //printf("endlinvra\n\n");
            //if(eq_orderreg[eqindx]>=0)break;
          }
          matsize=0;
          for (i=0; i<fdim; i++)if(dimmat[i]>matsize)matsize=dimmat[i];
          for (i=0; i<fdim; i++)if(dimmat[i]==matsize)break;
          //printf("eq %d i %d\n",eqindx,i);
          if(matsize==0) {
            eq_orderreg[eqindx]=tempint;
            ha_eqint[eqindx]=true;
          }
          else {
            eq_orderreg[eqindx]=i;
          }

          /*for (i=0; i<fdim; i++) {
            for(i4=0; i4<nlinvars; i4++) {
              for(i3=0; i3<ha_var[LinVars[i4].LinVarIndx].size; i3++) {
                if(strcmp(LinVars[i4].dimnames[i3],arSet[i].arIndx)==0&&ha_set[ha_eq[eqindx].setid[i]].regional&&orderreg[LinVars[i4].LinVarIndx]==i3) { //ha_set[ha_var[LinVars[i4].LinVarIndx].setid[i3]].regional) {
                  eq_orderreg[eqindx]=i;
                  break;
                }
              }
              if(eq_orderreg[eqindx]>=0)break;
            }
            if(eq_orderreg[eqindx]>=0)break;
          }
          if(eq_orderreg[eqindx]==-1){eq_orderreg[eqindx]=tempint;ha_eqint[eqindx]=true;}*/
          //printf("eq %s eqreg %d\n",ha_eq[eqindx].cofname,eq_orderreg[eqindx]);
        }

        /*for(i4=0; i4<nlinvars; i4++) {
          for(i3=0; i3<ha_var[LinVars[i4].LinVarIndx].size; i3++) {
            if(strcmp(LinVars[i4].dimnames[i3],arSet[eq_orderreg[eqindx]].arIndx)==0) {
              if(orderreg[LinVars[i4].LinVarIndx]==i3)varchoose[LinVars[i4].LinVarIndx]=true;
              else if(!varchoose[LinVars[i4].LinVarIndx]){
                printf("line %s\n",linecopy);
                printf("var %s varinx %d regindx %d vachoose %d eq indx %d i3 %d\n",ha_var[LinVars[i4].LinVarIndx].cofname,LinVars[i4].LinVarIndx,varchoose[LinVars[i4].LinVarIndx],orderreg[LinVars[i4].LinVarIndx],eq_orderreg[eqindx],i3);
                //orderreg[LinVars[i4].LinVarIndx]=i3;
                varchoose[LinVars[i4].LinVarIndx]=true;
              }
              break;
            }
          }
        }*/

        if(eq_orderreg[eqindx]>-1) {
          matsize=0;
          memset(linvarrcount,0,nlinvars*sizeof(bool));
          for(i4=0; i4<nlinvars; i4++) {
            for(lvar4=0; lvar4<i4; lvar4++)if(strcmp(LinVars[i4].LinVarName,LinVars[lvar4].LinVarName)==0)break;
            if(lvar4<i4&&linvarrcount[lvar4]) {
              linvarrcount[i4]=true;
              //printf("Here!!!!\n");
              continue;
            }
            //for(i3=0; i3<ha_var[LinVars[i4].LinVarIndx].size; i3++) {
            if(strcmp(LinVars[i4].dimnames[orderreg[LinVars[i4].LinVarIndx]],arSet[eq_orderreg[eqindx]].arIndx)!=0&&!var_inter[LinVars[i4].LinVarIndx]) {
              matsize+=ha_var[LinVars[i4].LinVarIndx].matsize;
              l01=0;
              for(j01=0; j01<ha_var[LinVars[i4].LinVarIndx].matsize; j01++)if(ha_cgeshock[ha_var[LinVars[i4].LinVarIndx].begadd+j01].ShockId)l01++;
              matsize+=ha_var[LinVars[i4].LinVarIndx].matsize-l01;
              linvarrcount[i4]=true;
              if(lvar4<i4)linvarrcount[lvar4]=true;
            }
            //}
          }
          if(matsize>=ha_eq[eqindx].matsize)ha_eqint[eqindx]=true;
          else {
            for(i4=0; i4<nlinvars; i4++)
              //for(i3=0; i3<ha_var[LinVars[i4].LinVarIndx].size; i3++) {
              if(strcmp(LinVars[i4].dimnames[orderreg[LinVars[i4].LinVarIndx]],arSet[eq_orderreg[eqindx]].arIndx)!=0&&!var_inter[LinVars[i4].LinVarIndx]) {
                var_inter[LinVars[i4].LinVarIndx]=true;
              }
            //}
          }
        }
        //printf("eq %s eqreg %d j %d arset %s\n",ha_eq[eqindx].cofname,eq_orderreg[eqindx],j,arSet[eq_orderreg[eqindx]].arIndx);
        /*if(j>0){
        //l=0;
        for(i4=0; i4<nlinvars; i4++){
          for(i3=0; i3<ha_var[LinVars[i4].LinVarIndx].size;i3++){
          if(strcmp(LinVars[i4].dimnames[i3],arSet[eq_orderreg[eqindx]].arIndx)==0){
            //printf("eqindx %d\n",eqindx);
            //printf("lanvar %s dimname %s i3 %d\n",ha_var[LinVars[i4].LinVarIndx].cofname,LinVars[i4].dimnames[i3],i3);
            if(orderreg[LinVars[i4].LinVarIndx]!=i3){var_inter[LinVars[i4].LinVarIndx]=true;break;}
            else break;
          }
          }
          if(i3==ha_var[LinVars[i4].LinVarIndx].size)var_inter[LinVars[i4].LinVarIndx]=true;
        }
        }*/
      }
      free(arSet);
      free(linvarrcount);
      //free(linvartcount);
      free(LinVars);
      eqindx++;
    }
  }
  //free(varchoose);
  fclose(filehandle);
  return 1;
}

int NestedMatvarRead(char *fname, char *commsyntax,ha_cgeset *ha_set,uvdim nset,ha_cgesetele *ha_setele,hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar,ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele,ha_cgeexovar *ha_cgeshock,bool *var_inter,hcge_cof *ha_eq,bool *ha_eqint,uvdim *eq_orderintra,uvdim *eq_orderreg,uvadd allregset,uvadd alltimeset,uvdim *orderintra,uvdim *orderreg) {
  FILE * filehandle;
  char tline[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[TABREADLINE],lintmp[TABREADLINE];//,*p1=NULL;
  char *readitem=NULL,*p=NULL,*p1=NULL;//,*p2=NULL,*varpnts;
  ha_cgetype zerodivide=0;
  uvdim fdim,np;
  uvadd rowindx,j,j01,l,l01,lj,dcountdim1[4*MAXVARDIM],dcountdim2[4*MAXVARDIM],nloops,nloopslin,li3,l2,matrow,eqindx=0,matsize;//,sizelinvars,totlinvars,templinvars
  uvadd sj,i,i3,i4,nlinvars,lvar,lvar1,lvar2,lvar3,lvar4,dcount,fdimlin=0,varindx1,varindx2,tempint;
  int leadlag;
  ha_cgetype dimmat[4*MAXVARDIM];
  //bool *varchoose= (bool *) calloc (nvar,sizeof(bool));
  //printf("varchoose %d\n",varchoose[70]);

  filehandle = fopen(fname,"r");
  matrow=0;
  //printf("OK!!!!\n");

  while (ha_cgertabl1(commsyntax,filehandle,line,ha_cofvar,ha_cof,ncof,&zerodivide,TABREADLINE)) {
    //PetscPrintf(PETSC_COMM_WORLD,"line %s\n",line);
    if (strstr(line,"(default")==NULL) {
      ha_cgefrstr1(line, commsyntax, "");
      ha_cgefrstr1(line, "(linear)", "");
      while (ha_cgefrstr(line,"  ", " "));
      while (ha_cgefrchr(line, '[', '('));
      while (ha_cgefrchr(line, ']', ')'));
      while (ha_cgefrchr(line, '{', '('));
      while (ha_cgefrchr(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=ha_cgenfind(line, "(all,");
      //printf("line %s\n",line);
      if (fdim==0) {
        readitem = strtok(line+1," ");
        strcpy(ha_eq[eqindx].cofname,readitem);
        //printf("eq name %s\n",readitem);
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      else {
        //printf("line %s\n",line+1);
        readitem = strtok(line+1,"(");
        //printf("read %s\n",readitem);
        strcpy(ha_eq[eqindx].cofname,readitem);
        //printf("eq name %s\n",readitem);
        strcpy(line,linecopy);
        i=ha_cgerevfind(line, "(all,");
        readitem=line+i;
        readitem = strtok(readitem,")");
        readitem = strtok(NULL,"=");
        //printf("left %s\n",readitem);
        strcpy(vname,readitem);
        //printf("left %s linecopy %s\n",linecopy1,linecopy);
        strcpy(line,linecopy);
        //printf("left %s\n",line);
        readitem = strtok(line,"=");
        //printf("read %s\n",readitem);
        readitem = strtok(NULL,";");
        //printf("left %s\n",readitem);
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      while (ha_cgefrstr(readitem," ", ""));
      while (ha_cgerecovar(readitem)==1);
      hnew_intrpl(readitem);
      //printf("read %s\n",readitem);
      strcpy(tline,readitem);
      strcpy(line1,readitem);
      np=ha_cgenfind(readitem,"p_");
      hcge_linvars *LinVars= (hcge_linvars *) calloc (np+1,sizeof(hcge_linvars));
      i3=0;
      //sizelinvars=0;
      lvar=0;
      for (i=0; i<np; i++) {
        //varindx1=0;
        varindx2=0;
        while(-1<0) {
          varindx1=ha_cgefind(readitem+varindx2,"p_");
          if(varindx1==-1) break;
          varindx2=varindx2+varindx1;
          if(varindx2==0||readitem[varindx2-1]=='*'||readitem[varindx2-1]=='+'||readitem[varindx2-1]=='-'||readitem[varindx2-1]=='('||readitem[varindx2-1]==',') break;
          else varindx2++;
        }
        if(varindx1==-1) break;
        lvar=lvar+varindx2;
        readitem=readitem+varindx2+2;
        //printf("varbeg %s\n",readitem);
        p=strpbrk(readitem,"{+*-/^)");
        if(p!=NULL) {
          strncpy(vname,readitem,p-readitem);
          vname[p-readitem]='\0';
          //printf("vname %s\n",vname);
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(ha_var[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          p=strpbrk(readitem,"}+*-/^)");
          if (*p=='}') {
            p++;
          }
          strncpy(tline,readitem,p-readitem);
          tline[p-readitem]='\0';
          switch (ha_var[l].size) {
          case 0:
            break;
          case 1:
            p = strtok(tline,"{");
            p = strtok(NULL,"}");
            leadlag=0;
            hnew_arset(p,&leadlag);
            if(leadlag!=0) {
              var_inter[l]=true;
              //printf("var!!!!!!!!!!!!!!!!!!! %s\n",ha_var[l].cofname);
              orderintra[l]=-1;
              orderreg[l]=-1;
            }//printf("var %s\n",ha_var[l].cofname);}
            strcpy(LinVars[i3].dimnames[0],p);
            //LinVars[i3].dimleadlag[0]=leadlag;
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l01=ha_cgefind(linecopy,lintmp);
            //printf("lin %s\n",lintmp);
            if (l01>-1) {
              p1=&linecopy[0]+l01;
              strncpy(LinVars[i3].dimsetnames[0],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              //printf("lin1 %s\n",lintmp);
              lvar1=ha_cgenfind(linecopy,lintmp);
              lvar3=ha_cgefind(linecopy,lintmp);
              //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
                  if (lvar4>-1&&lvar4<lvar) {
                    lvar3=lvar3+lvar4+4;
                  }
                  else {
                    break;
                  }
                }
              //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
              p1=&linecopy[0]+lvar3;
              //p1=&linecopy[0]+ha_cgefind(linecopy,lintmp);
              p1=p1+strlen(lintmp);
              //printf("lin %s\n",lintmp);
              strncpy(LinVars[i3].dimsetnames[0],p1,strchr(p1,',')-p1);
            }
            //printf("vname %s lin %s\n",vname,LinVars[i3].dimsetnames[0]);
            //printf("lin %s indx %s\n",ha_var[l].cofname,p);
            break;
          default:
            p = strtok(tline,"{");
            for (i4=0; i4<ha_var[l].size-1; i4++) {
              p = strtok(NULL,",");
              leadlag=0;
              hnew_arset(p,&leadlag);
              if(leadlag!=0) {
                var_inter[l]=true;
                //printf("var!!!!!!!!!!!!!!!!!!! %s\n",ha_var[l].cofname);
                orderintra[l]=-1;
                orderreg[l]=-1;
              }//printf("var %s\n",ha_var[l].cofname);}
              strcpy(LinVars[i3].dimnames[i4],p);
              //printf("p %s\n",p);
              //LinVars[i3].dimleadlag[i4]=leadlag;
              strcpy(lintmp,"(all,");
              strcat(lintmp,p);
              strcat(lintmp,",");
              l01=ha_cgefind(linecopy,lintmp);
              //printf("lin %s\n",lintmp);
              if (l01>-1) {
                p1=&linecopy[0]+l01;
                strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
              }
              else {
                strcpy(lintmp,"sum(");
                strcat(lintmp,p);
                strcat(lintmp,",");
                //printf("lin1 %s\n",lintmp);
                lvar1=ha_cgenfind(linecopy,lintmp);
                lvar3=ha_cgefind(linecopy,lintmp);
                //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
                if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                    lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
                    if (lvar4>-1&&lvar4<lvar) {
                      lvar3=lvar3+lvar4+4;
                    }
                    else {
                      break;
                    }
                  }
                //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
                p1=&linecopy[0]+lvar3;
                //p1=&linecopy[0]+ha_cgefind(linecopy,lintmp);
                p1=p1+strlen(lintmp);
                //printf("lin %s p1 %s len %d\n",lintmp,p1,strchr(p1,',')-p1-strlen(lintmp));
                strncpy(LinVars[i3].dimsetnames[i4],p1,strchr(p1,',')-p1);
                //printf("lin %s\n",LinVars[i3].dimsetnames[i4]);
              }
              //printf("vname %s lin %s\n",vname,LinVars[i3].dimsetnames[i4]);
            }
            p = strtok(NULL,"}");
            leadlag=0;
            hnew_arset(p,&leadlag);
            if(leadlag!=0) {
              var_inter[l]=true;
              //printf("var!!!!!!!!!!!!!!!!!!! %s\n",ha_var[l].cofname);
              orderintra[l]=-1;
              orderreg[l]=-1;
            }//printf("var %s\n",ha_var[l].cofname);}
            strcpy(LinVars[i3].dimnames[i4],p);
            //LinVars[i3].dimleadlag[i4]=leadlag;
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l01=ha_cgefind(linecopy,lintmp);
            //printf("lin %s\n",lintmp);
            if (l01>-1) {
              p1=&linecopy[0]+l01;
              strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=ha_cgenfind(linecopy,lintmp);
              lvar3=ha_cgefind(linecopy,lintmp);
              //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
                  if (lvar4>-1&&lvar4<lvar) {
                    lvar3=lvar3+lvar4+4;
                  }
                  else {
                    break;
                  }
                }
              //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
              p1=&linecopy[0]+lvar3;
              p1=p1+strlen(lintmp);
              //printf("lin %s p1 %s\n",lintmp,p1);
              strncpy(LinVars[i3].dimsetnames[i4],p1,strchr(p1,',')-p1);
            }
            //printf("vname %s lin %s l1 %d\n",vname,LinVars[i3].dimsetnames[i4],l1);
            break;
          }
          i3++;
          //}
          //LinVars[i][p-readitem]='\0';
        }
        else {
          strcpy(vname,readitem);
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(ha_var[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          i3++;
          //}
          //break;
        }
      }
      nlinvars=i3;
      for (i=0; i<nlinvars; i++) i3+=ha_var[LinVars[i].LinVarIndx].size;
      ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (fdim+i3+1,sizeof(ha_cgesetindx));
      bool *linvarrcount= (bool *) calloc (nlinvars,sizeof(bool));

      strcpy(line,linecopy);
      //printf("line %s\n",line);
      nloops=1;
      if (fdim>0) {
        for (i=0; i<fdim; i++) {
          if(i==0) {
            readitem = strtok(line,",");
          }
          else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
          strcpy(arSet[i].arIndx,readitem);
          readitem = strtok(NULL,")");
          for (i4=0; i4<nset; i4++) if(strcmp(readitem,ha_set[i4].setname)==0) {
              arSet[i].setid=i4;
              ha_eq[eqindx].setid[i]=i4;
              break;
            }
          nloops=nloops*ha_set[arSet[i].setid].size;
        }
      }
      ha_eq[eqindx].size=fdim;
      ha_eq[eqindx].matsize=nloops;


      if(eqindx==0) {
        for(i=0; i<nvar; i++) {
          if(var_inter[i]) continue;
          j=0;
          for(i4=0; i4<ha_var[i].size; i4++)if(ha_set[ha_var[i].setid[i4]].intertemp) {
              orderintra[i]=i4;
              j++;
              break;
            }
          if(j==0) {
            var_inter[i]=true;
            //printf("var!!!!!!!!!!!!!!!!!!! %s\n",ha_var[i].cofname);
            orderreg[i]=-1;
            orderintra[i]=-1;
          }
          else {
            j01=0;
            for(i4=0; i4<ha_var[i].size; i4++)if(ha_set[ha_var[i].setid[i4]].regional) {
                orderreg[i]=i4;
                j01++;
                //break;
              }
            //if(j1>1)
            //printf("var %s j01 %d order %d\n",ha_var[i].cofname,j01,orderreg[i]);
            //if(i==23)orderreg[i]=1;
          }
        }
      }


      j=0;
      for(i4=0; i4<nlinvars; i4++)for(i=0; i<ha_var[LinVars[i4].LinVarIndx].size; i++) if(ha_set[ha_var[LinVars[i4].LinVarIndx].setid[i]].intertemp)j++;
      if(j==0) {
        ha_eqint[eqindx]=true;//eq_orderreg[eqindx]=-1;eq_orderintra[eqindx]=-1;}
        //printf("Equation %ld\n",eqindx);
      }
      j=0;
      for (i=0; i<fdim; i++)if(ha_set[ha_eq[eqindx].setid[i]].regional) {
          j++;
          eq_orderreg[eqindx]=i;
        }
      if(j>=2) {
        //tempint=eq_orderreg[eqindx];
        //eq_orderreg[eqindx]=-1;
        for (i=0; i<fdim; i++)dimmat[i]=-1;
        for (i=0; i<fdim; i++) {
          memset(linvarrcount,0,nlinvars*sizeof(bool));
          for(i4=0; i4<nlinvars; i4++) {
            for(lvar4=0; lvar4<i4; lvar4++)if(strcmp(LinVars[i4].LinVarName,LinVars[lvar4].LinVarName)==0)break;
            if(lvar4<i4&&linvarrcount[lvar4]) {
              linvarrcount[i4]=true;
              continue;
            }
            for(i3=0; i3<ha_var[LinVars[i4].LinVarIndx].size; i3++) {
              if(strcmp(LinVars[i4].dimnames[i3],arSet[i].arIndx)==0&&ha_set[ha_eq[eqindx].setid[i]].regional&&orderreg[LinVars[i4].LinVarIndx]==i3) {
                dimmat[i]+=ha_var[LinVars[i4].LinVarIndx].matsize;
                linvarrcount[i4]=true;
                if(lvar4<i4)linvarrcount[lvar4]=true;
                l01=0;
                for(j01=0; j01<ha_var[LinVars[i4].LinVarIndx].matsize; j01++)if(ha_cgeshock[ha_var[LinVars[i4].LinVarIndx].begadd+j01].ShockId)l01++;
                dimmat[i]-=l01;
              }
            }
          }
        }
        matsize=0;
        for (i=0; i<fdim; i++)if(dimmat[i]>matsize)matsize=dimmat[i];
        for (i=0; i<fdim; i++)if(dimmat[i]==matsize)break;
        if(matsize==0) {
          eq_orderreg[eqindx]=-1;
        }
        else {
          eq_orderreg[eqindx]=i;
        }
      }
      if(eq_orderreg[eqindx]>-1) {
        matsize=0;
        memset(linvarrcount,0,nlinvars*sizeof(bool));
        for(i4=0; i4<nlinvars; i4++) {
          for(lvar4=0; lvar4<i4; lvar4++)if(strcmp(LinVars[i4].LinVarName,LinVars[lvar4].LinVarName)==0)break;
          if(lvar4<i4&&linvarrcount[lvar4]) {
            linvarrcount[i4]=true;
            continue;
          }
          if(strcmp(LinVars[i4].dimnames[orderreg[LinVars[i4].LinVarIndx]],arSet[eq_orderreg[eqindx]].arIndx)!=0&&orderreg[LinVars[i4].LinVarIndx]!=-1) {
//              j=0;
//              for(i=0; i<ha_var[LinVars[i4].LinVarIndx].size; i++)if(strcmp(LinVars[i4].dimnames[i],arSet[eq_orderreg[eqindx]].arIndx)==0)j++;
//              if(j==0)
            l01=0;
            for(j01=0; j01<ha_var[LinVars[i4].LinVarIndx].matsize; j01++)if(ha_cgeshock[ha_var[LinVars[i4].LinVarIndx].begadd+j01].ShockId)l01++;
            matsize+=ha_var[LinVars[i4].LinVarIndx].matsize-l01;
            linvarrcount[i4]=true;
            if(lvar4<i4)linvarrcount[lvar4]=true;
          }
        }
        if(matsize>=ha_eq[eqindx].matsize)eq_orderreg[eqindx]=-1;//printf("OOOOOOOOOOOOOOOOO\n");}
        else {
          for(i4=0; i4<nlinvars; i4++)
            if(strcmp(LinVars[i4].dimnames[orderreg[LinVars[i4].LinVarIndx]],arSet[eq_orderreg[eqindx]].arIndx)!=0&&orderreg[LinVars[i4].LinVarIndx]!=-1) {
              //var_inter[LinVars[i4].LinVarIndx]=true;
//              j=0;
//              for(i=0; i<ha_var[LinVars[i4].LinVarIndx].size; i++)if(strcmp(LinVars[i4].dimnames[i],arSet[eq_orderreg[eqindx]].arIndx)==0)j++;
//              if(j==0)
              orderreg[LinVars[i4].LinVarIndx]=-1;//orderintra[LinVars[i4].LinVarIndx]=-1;
            }
        }
      }


      j=0;
      for (i=0; i<fdim; i++)if(ha_set[ha_eq[eqindx].setid[i]].intertemp) {
          j++;
          eq_orderintra[eqindx]=i;
        }
      if(j==0) {
        //for (i4=0; i4<nlinvars; i4++) var_inter[LinVars[i4].LinVarIndx]=true;
        ha_eqint[eqindx]=true;
        //printf("Equation1 %ld\n",eqindx);
        //eq_orderreg[eqindx]=-1;
        //eq_orderintra[eqindx]=-1;
      }

      if(ha_eqint[eqindx]==true) {
        eq_orderreg[eqindx]=-1;
        eq_orderintra[eqindx]=-1;
      }

      free(arSet);
      free(LinVars);
      free(linvarrcount);
      eqindx++;
    }
  }
  //free(varchoose);
  fclose(filehandle);
  return 1;
}

int NewMatreadele(char *fname, char *commsyntax,ha_cgeset *ha_set,uvdim nset,ha_cgesetele *ha_setele,hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar,ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele, uvadd nexo,ha_cgeexovar *ha_cgeshock,uvadd ndblock,uvadd alltimeset,uvadd allregset,bool *ha_eqint,uvadd *ha_eqadd,uvdim *ha_eqtime,uvdim *ha_eqreg,uvadd *counteq,uvadd nintraeq,bool *sbbd_overrid,PetscInt Istart,PetscInt Iend,PetscInt *dnz,PetscInt *dnnz,PetscInt *onz,PetscInt *onnz,PetscInt *dnzB,PetscInt *dnnzB,PetscInt *onzB,PetscInt *onnzB,int nesteddbbd) {
  FILE * filehandle;
  char tline[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[TABREADLINE],lintmp[TABREADLINE];//,*p1=NULL;
  char *readitem=NULL,*p=NULL,*p1=NULL;//,*p2=NULL,*varpnts;
  PetscInt Iindx=0,Jindx;
  ha_cgetype zerodivide=0;
  uvdim fdim=0,np,i4,sup,supset[MAXSUPSET];
  uvadd rowindx,rowindxorg,l,l1,lj,dcountdim1[4*MAXVARDIM],dcountdim2[4*MAXVARDIM],dcountdim3[4*MAXVARDIM],dcountdim4[4*MAXVARDIM],dcountdim5[4*MAXVARDIM],nloops,nloopslin,nloopsfac,li3,l2,matrow,matroworg,ltime,lreg,leq=0,eqindx=0;//,sizelinvars,totlinvars,templinvars
  uvadd nreg=0,nint=0,sj,i,i3;
  if(allregset>-1)nreg=ha_set[allregset].size;
  if(alltimeset>-1)nint=ha_set[alltimeset].size;
  int nlinvars,lvar,lvar1,lvar2,lvar3,lvar4,dcount,fdimlin=0,leadlag,varindx1,varindx2;
  uvadd *counteq1= (uvadd *) calloc (ndblock,sizeof(uvadd));
  for(i=0; i<ndblock; i++)counteq1[i]=counteq[i];
  filehandle = fopen(fname,"r");
  matrow=nintraeq;
  matroworg=0;
  Jindx=0;
  //printf("OK!!!!\n");

  while (ha_cgertabl1(commsyntax,filehandle,line,ha_cofvar,ha_cof,ncof,&zerodivide,TABREADLINE)) {
    //PetscPrintf(PETSC_COMM_WORLD,"line %s\n",line);
    if (strstr(line,"(default")==NULL) {
      ha_cgefrstr1(line, commsyntax, "");
      ha_cgefrstr1(line, "(linear)", "");
      while (ha_cgefrstr(line,"  ", " "));
      while (ha_cgefrchr(line, '[', '('));
      while (ha_cgefrchr(line, ']', ')'));
      while (ha_cgefrchr(line, '{', '('));
      while (ha_cgefrchr(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=ha_cgenfind(line, "(all,");
      //printf("line %s\n",line);
      if (fdim==0) {
        readitem = strtok(line+1," ");
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      else {
        i=ha_cgerevfind(line, "(all,");
        readitem=line+i;
        readitem = strtok(readitem,")");
        readitem = strtok(NULL,"=");
        //printf("left %s\n",readitem);
        strcpy(vname,readitem);
        //printf("left %s linecopy %s\n",linecopy1,linecopy);
        strcpy(line,linecopy);
        //printf("left %s\n",line);
        readitem = strtok(line,"=");
        //printf("read %s\n",readitem);
        readitem = strtok(NULL,";");
        //printf("left %s\n",readitem);
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      while (ha_cgefrstr(readitem," ", ""));
      while (ha_cgerecovar(readitem)==1);
      hnew_intrpl(readitem);
      //printf("read %s\n",readitem);
      strcpy(tline,readitem);
      strcpy(line1,readitem);
      np=ha_cgenfind(readitem,"p_");
      hcge_linvars *LinVars= (hcge_linvars *) calloc (np+1,sizeof(hcge_linvars));
      i3=0;
      //sizelinvars=0;
      lvar=0;
      for (i=0; i<np; i++) {
        //varindx1=0;
        varindx2=0;
        while(-1<0) {
          varindx1=ha_cgefind(readitem+varindx2,"p_");
          if(varindx1==-1) break;
          varindx2=varindx2+varindx1;
          if(varindx2==0||readitem[varindx2-1]=='*'||readitem[varindx2-1]=='+'||readitem[varindx2-1]=='-'||readitem[varindx2-1]=='('||readitem[varindx2-1]==',') break;
          else varindx2++;
        }
        if(varindx1==-1) break;
        lvar=lvar+varindx2;
        readitem=readitem+varindx2+2;
        //printf("varbeg %s\n",readitem);
        p=strpbrk(readitem,"{+*-/^)");
        if(p!=NULL) {
          strncpy(vname,readitem,p-readitem);
          vname[p-readitem]='\0';
          //printf("vname %s\n",vname);
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(ha_var[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          p=strpbrk(readitem,"}+*-/^)");
          if (*p=='}') {
            p++;
          }
          strncpy(tline,readitem,p-readitem);
          tline[p-readitem]='\0';
          switch (ha_var[l].size) {
          case 0:
            break;
          case 1:
            p = strtok(tline,"{");
            p = strtok(NULL,"}");
            leadlag=0;
            hnew_arset(p,&leadlag);
            strcpy(LinVars[i3].dimnames[0],p);
            LinVars[i3].dimleadlag[0]=leadlag;
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=ha_cgefind(linecopy,lintmp);
            //printf("lin %s\n",lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[0],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              //printf("lin1 %s\n",lintmp);
              lvar1=ha_cgenfind(linecopy,lintmp);
              lvar3=ha_cgefind(linecopy,lintmp);
              //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
                  if (lvar4>-1&&lvar4<lvar) {
                    lvar3=lvar3+lvar4+4;
                  }
                  else {
                    break;
                  }
                }
              //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
              p1=&linecopy[0]+lvar3;
              //p1=&linecopy[0]+ha_cgefind(linecopy,lintmp);
              p1=p1+strlen(lintmp);
              //printf("lin %s\n",lintmp);
              strncpy(LinVars[i3].dimsetnames[0],p1,strchr(p1,',')-p1);
            }
            //printf("vname %s lin %s\n",vname,LinVars[i3].dimsetnames[0]);
            //printf("lin %s indx %s\n",ha_var[l].cofname,p);
            break;
          default:
            p = strtok(tline,"{");
            for (i4=0; i4<ha_var[l].size-1; i4++) {
              p = strtok(NULL,",");
              leadlag=0;
              hnew_arset(p,&leadlag);
              strcpy(LinVars[i3].dimnames[i4],p);
              LinVars[i3].dimleadlag[i4]=leadlag;
              strcpy(lintmp,"(all,");
              strcat(lintmp,p);
              strcat(lintmp,",");
              l1=ha_cgefind(linecopy,lintmp);
              //printf("lin %s\n",lintmp);
              if (l1>-1) {
                p1=&linecopy[0]+l1;
                strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
              }
              else {
                strcpy(lintmp,"sum(");
                strcat(lintmp,p);
                strcat(lintmp,",");
                //printf("lin1 %s\n",lintmp);
                lvar1=ha_cgenfind(linecopy,lintmp);
                lvar3=ha_cgefind(linecopy,lintmp);
                //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
                if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                    lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
                    if (lvar4>-1&&lvar4<lvar) {
                      lvar3=lvar3+lvar4+4;
                    }
                    else {
                      break;
                    }
                  }
                //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
                p1=&linecopy[0]+lvar3;
                //p1=&linecopy[0]+ha_cgefind(linecopy,lintmp);
                p1=p1+strlen(lintmp);
                //printf("lin %s p1 %s len %d\n",lintmp,p1,strchr(p1,',')-p1-strlen(lintmp));
                strncpy(LinVars[i3].dimsetnames[i4],p1,strchr(p1,',')-p1);
                //printf("lin %s\n",LinVars[i3].dimsetnames[i4]);
              }
              //printf("vname %s lin %s\n",vname,LinVars[i3].dimsetnames[i4]);
            }
            p = strtok(NULL,"}");
            leadlag=0;
            hnew_arset(p,&leadlag);
            strcpy(LinVars[i3].dimnames[i4],p);
            LinVars[i3].dimleadlag[i4]=leadlag;
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=ha_cgefind(linecopy,lintmp);
            //printf("lin %s\n",lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=ha_cgenfind(linecopy,lintmp);
              lvar3=ha_cgefind(linecopy,lintmp);
              //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
                  if (lvar4>-1&&lvar4<lvar) {
                    lvar3=lvar3+lvar4+4;
                  }
                  else {
                    break;
                  }
                }
              //printf("lvar1 %d linecopy %s\n",lvar1,linecopy);
              p1=&linecopy[0]+lvar3;
              p1=p1+strlen(lintmp);
              //printf("lin %s p1 %s\n",lintmp,p1);
              strncpy(LinVars[i3].dimsetnames[i4],p1,strchr(p1,',')-p1);
            }
            //printf("vname %s lin %s l1 %d\n",vname,LinVars[i3].dimsetnames[i4],l1);
            break;
          }
          i3++;
          //}
          //LinVars[i][p-readitem]='\0';
        }
        else {
          strcpy(vname,readitem);
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(ha_var[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          i3++;
          //}
          //break;
        }
      }
      nlinvars=i3;
      for (i=0; i<nlinvars; i++) i3+=ha_var[LinVars[i].LinVarIndx].size;
      ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (fdim+i3+1,sizeof(ha_cgesetindx));

      strcpy(line,linecopy);
      //printf("line12 %s\n",line);
      nloops=1;
      if (fdim>0) {
        for (i=0; i<fdim; i++) {
          if(i==0) {
            readitem = strtok(line,",");
          }
          else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
          strcpy(arSet[i].arIndx,readitem);
          readitem = strtok(NULL,")");
          //strcpy(arSet[i].arSet,readitem);
          for (i4=0; i4<nset; i4++) if(strcmp(readitem,ha_set[i4].setname)==0) {
              arSet[i].setid=i4;
              if(ha_set[i4].intertemp) *sbbd_overrid=true;
              //arSet[i].SetSize=ha_set[i4].size;
              //arSet[i].SetBegAdd=ha_set[i4].begadd;
              //arSet[i].subsetid=ha_set[i4].subsetid;
              //if(ha_set[i4].subsetid==1) {
              //arSet[i].SuperSetSize=ha_set[i4].supersetsize;
              //check=true;
              //} else {
              //arSet[i].SuperSetSize=ha_set[i4].size;
              //}
              break;
            }
          nloops=nloops*ha_set[arSet[i].setid].size;
        }
        //*****Adjust ha_eqreg*******************//
        /*if(allregset>=0)for (i=0; i<nlinvars; i++)for (i4=0; i4<ha_var[LinVars[i].LinVarIndx].size; i4++)
          if(ha_var[LinVars[i].LinVarIndx].setid[i4]==allregset)
            for (dcount=0; dcount<fdim; dcount++)
              if(strcmp(LinVars[i].dimnames[i4],arSet[dcount].arIndx)==0){printf("line %slinvar %s arset %d varset %d\n\n",linecopy,ha_var[LinVars[i].LinVarIndx].cofname,dcount,i4);ha_eqreg[eqindx]=dcount;break;}*/
        //*****End Adjust ha_eqreg*******************//
        dcountdim1[fdim-1]=1;
        for (dcount=fdim-2; dcount>-1; dcount--) {
          dcountdim1[dcount]=ha_set[arSet[dcount+1].setid].size*dcountdim1[dcount+1];
        }
        if(alltimeset>=0&&allregset<0) {
          for (dcount=ha_eqtime[eqindx]-1; dcount>-1; dcount--) {
            dcountdim3[dcount]=dcountdim1[dcount]/ha_set[arSet[ha_eqtime[eqindx]].setid].size;
          }
          for (dcount=fdim-1; dcount>ha_eqtime[eqindx]-1; dcount--) {
            dcountdim3[dcount]=dcountdim1[dcount];
          }
        }
        if(alltimeset<0&&allregset>=0) {
          for (dcount=ha_eqreg[eqindx]-1; dcount>-1; dcount--) {
            dcountdim4[dcount]=dcountdim1[dcount]/ha_set[arSet[ha_eqreg[eqindx]].setid].size;
          }
          for (dcount=fdim-1; dcount>ha_eqreg[eqindx]-1; dcount--) {
            dcountdim4[dcount]=dcountdim1[dcount];
          }
        }
        if(alltimeset>=0&&allregset>=0) {
          i3=1;
          dcountdim4[fdim-1]=dcountdim1[fdim-1];
          for (dcount=fdim-2; dcount>-1; dcount--) {
            if(dcount+1==ha_eqreg[eqindx])i3*=ha_set[arSet[ha_eqreg[eqindx]].setid].size;
            if(dcount+1==ha_eqtime[eqindx])i3*=ha_set[arSet[ha_eqtime[eqindx]].setid].size;
            dcountdim4[dcount]=dcountdim1[dcount]/i3;
          }
        }
      }
      //printf("OK!!!\n");
      for (i=0; i<nlinvars; i++) {
        i3=0;
        nloopslin=nloops;
        if (fdim==0) {
          for (l2=0; l2<ha_var[LinVars[i].LinVarIndx].size; l2++) {
            strcpy(arSet[i3].arIndx,LinVars[i].dimnames[l2]);
            //strcpy(arSet[i3].arSet,LinVars[i].dimsetnames[l2]);
            for (i4=0; i4<nset; i4++) if(strcmp(ha_set[i4].setname,LinVars[i].dimsetnames[l2])==0) {
                break;
              }
            arSet[i3].setid=i4;
            //arSet[i3].SetSize=ha_set[i4].size;
            //arSet[i3].SetBegAdd=ha_set[i4].begadd;
            nloopslin=nloopslin*ha_set[arSet[i3].setid].size;
            if(ha_set[i4].intertemp)*sbbd_overrid=true;
            //printf("set %s nloop %d nlooplin %d *sbbd_overrid %d\n",ha_set[i4].setname,nloops,nloopslin,*sbbd_overrid);
            i3++;
          }
        }
        else for (l2=0; l2<ha_var[LinVars[i].LinVarIndx].size; l2++) {
            for (sj=0; sj<fdim; sj++) {
              if(strcmp(arSet[sj].arIndx,LinVars[i].dimnames[l2])==0) {
                break;
              }
              else if (sj==fdim-1) {
                strcpy(arSet[fdim+i3].arIndx,LinVars[i].dimnames[l2]);
                //strcpy(arSet[fdim+i3].arSet,LinVars[i].dimsetnames[l2]);
                for (i4=0; i4<nset; i4++) if(strcmp(ha_set[i4].setname,LinVars[i].dimsetnames[l2])==0) {
                    break;
                  }
                arSet[fdim+i3].setid=i4;
                //arSet[fdim+i3].SetSize=ha_set[i4].size;
                //arSet[fdim+i3].SetBegAdd=ha_set[i4].begadd;
                nloopslin=nloopslin*ha_set[arSet[fdim+i3].setid].size;
                if(ha_set[i4].intertemp)*sbbd_overrid=true;
                //printf("set %s nloop %d nlooplin %d *sbbd_overrid %d\n",ha_set[i4].setname,nloops,nloopslin,*sbbd_overrid);
                i3++;
              }
            }
          }
        fdimlin=fdim+i3;
        //printf("userid %d\n",*sbbd_overrid);
        //printf("fdim %d fdimlin %d nloop %d nlooplin %d\n",fdim,fdimlin,nloops,nloopslin);
        //for (i4=0;i4<fdimlin;i4++) printf("indx1 %s set %s\n",arSet[i4].arIndx,arSet[i4].arSet);
        if (i3>0) {
          dcountdim2[fdimlin-1]=1;
          //printf("set %s\n",arSet[fdimlin-1].arSet);
          for (dcount=fdimlin-2; dcount>-1; dcount--) {
            //printf("set %s\n",arSet[dcount].arSet);
            dcountdim2[dcount]=ha_set[arSet[dcount+1].setid].size*dcountdim2[dcount+1];
          }
        }
        else for (dcount=0; dcount<fdim; dcount++) {
            dcountdim2[dcount]=dcountdim1[dcount];
          }
        nloopsfac=(uvadd)nloopslin/nloops;
        for (dcount=0; dcount<ha_var[LinVars[i].LinVarIndx].size; dcount++) {
          for (i4=0; i4<fdimlin; i4++) {
            if (strcmp(LinVars[i].dimnames[dcount],arSet[i4].arIndx)==0) {
              dcountdim5[dcount]=i4;
              break;
            }
          }
        }
        for(dcount=0; dcount<MAXSUPSET; dcount++)supset[dcount]=0;
        //printf("OKKKK3\n");
        for (dcount=0; dcount<ha_var[LinVars[i].LinVarIndx].size; dcount++) {
          if(ha_set[ha_var[LinVars[i].LinVarIndx].setid[dcount]].size!=ha_set[arSet[dcountdim5[dcount]].setid].size) {
            //printf("nrow %d vaset %s set ste %d\n",matrow,ha_set[ha_var[LinVars[i].LinVarIndx].setid[dcount]].setname,ha_set[arSet[dcountdim3[dcount]].setid].subsetid[1]);
            for(sup=1; sup<MAXSUPSET; sup++)if(ha_var[LinVars[i].LinVarIndx].setid[dcount]==ha_set[arSet[dcountdim5[dcount]].setid].subsetid[sup]) {
                supset[dcount]=sup;
                break;
              }
          }
          else supset[dcount]=0;
        }
        //printf("OKKKK3\n");
        for (lj=0; lj<nloopslin; lj++) {
          //if(lj==0)printf("OKKKK3z\n");
          l2=lj;
          rowindx=0;
          rowindxorg=0;
          for (dcount=0; dcount<fdimlin; dcount++) {
            l1=(uvadd) l2/dcountdim2[dcount];
            arSet[dcount].indx=l1;
            if(ha_eqint[eqindx]) {
              //if (dcount<fdim)rowindxorg=rowindxorg+l1*dcountdim1[dcount];
              if(alltimeset>=0&&allregset<0) {
                if (dcount<fdim&&dcount!=ha_eqtime[eqindx]) {
                  rowindx=rowindx+l1*dcountdim3[dcount];
                }
                if(dcount==ha_eqtime[eqindx])ltime=l1;
              }
              if(alltimeset<0&&allregset>=0) {
                if (dcount<fdim&&dcount!=ha_eqreg[eqindx]) {
                  rowindx=rowindx+l1*dcountdim4[dcount];
                }
                if(dcount==ha_eqreg[eqindx])ltime=l1;
              }
              if(alltimeset>=0&&allregset>=0) {
                if(nesteddbbd==1) {
                  if(ha_eqreg[eqindx]>-1) {
                    if (dcount<fdim&&dcount!=ha_eqtime[eqindx]&&dcount!=ha_eqreg[eqindx])rowindx=rowindx+l1*dcountdim4[dcount];
                  }
                  else {
                    if (dcount<fdim&&dcount!=ha_eqtime[eqindx])rowindx=rowindx+l1*dcountdim4[dcount];
                  }
                  if(dcount==ha_eqtime[eqindx])ltime=l1;
                  if(dcount==ha_eqreg[eqindx])lreg=l1;
                }
                else {
                  if (dcount<fdim&&dcount!=ha_eqtime[eqindx]&&dcount!=ha_eqreg[eqindx]) {
                    rowindx=rowindx+l1*dcountdim4[dcount];
                  }
                  if(dcount==ha_eqtime[eqindx])ltime=l1;
                  if(dcount==ha_eqreg[eqindx])lreg=l1;
                }
              }
            }
            else {
              if (dcount<fdim) {
                rowindx=rowindx+l1*dcountdim1[dcount];
                //rowindxorg=rowindxorg+l1*dcountdim1[dcount];
              }
            }
            //if (dcount<fdim) {
            //rowindx=rowindx+l1*dcountdim1[dcount];
            //}
            l2=l2-l1*dcountdim2[dcount];
          }
          li3=0;
          for (dcount=0; dcount<ha_var[LinVars[i].LinVarIndx].size; dcount++) {
            if(supset[dcount]==0) {
              li3=li3+(arSet[dcountdim5[dcount]].indx+LinVars[i].dimleadlag[dcount])*ha_var[LinVars[i].LinVarIndx].antidims[dcount];
            }
            else {
              li3=li3+(ha_setele[ha_set[arSet[dcountdim5[dcount]].setid].begadd+arSet[dcountdim5[dcount]].indx].setsh[supset[dcount]]+LinVars[i].dimleadlag[dcount])*ha_var[LinVars[i].LinVarIndx].antidims[dcount];
            }
          }

          /*for (dcount=0; dcount<ha_var[LinVars[i].LinVarIndx].size; dcount++) {
            for (i4=0; i4<fdimlin; i4++) {
              if (strcmp(LinVars[i].dimnames[dcount],arSet[i4].arIndx)==0) {
                if(ha_set[ha_var[LinVars[i].LinVarIndx].setid[dcount]].size==ha_set[arSet[i4].setid].size) {
                  li3=li3+(arSet[i4].indx+LinVars[i].dimleadlag[dcount])*ha_var[LinVars[i].LinVarIndx].antidims[dcount];
                } else {
                  for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[i4].setid].subsetid[sup]==ha_var[LinVars[i].LinVarIndx].setid[dcount]){
                    li3=li3+(ha_setele[ha_set[arSet[i4].setid].begadd+arSet[i4].indx].setsh[sup]+LinVars[i].dimleadlag[dcount])*ha_var[LinVars[i].LinVarIndx].antidims[dcount];
                    break;
                  }
                }
                break;
              }
            }
          }*/
          if(ha_eqint[eqindx]) {
            if(alltimeset>=0&&allregset>=0) {
              if(nesteddbbd==1) {
                if(ha_eqreg[eqindx]>-1)Jindx=counteq1[(ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+ltime].setsh[ha_set[arSet[ha_eqtime[eqindx]].setid].intsup])*(nreg+1)+ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+lreg].setsh[ha_set[arSet[ha_eqreg[eqindx]].setid].regsup]]+rowindx;
                else Jindx=counteq1[(ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+ltime].setsh[ha_set[arSet[ha_eqtime[eqindx]].setid].intsup])*(nreg+1)+nreg]+rowindx;
//if(Jindx==22439)printf("setname %s setele %d\n",ha_set[arSet[ha_eqtime[eqindx]].setid].setname,ha_set[arSet[ha_eqtime[eqindx]].setid].intsup);
              }
              else {
                if(arSet[ha_eqtime[eqindx]].setid==alltimeset)
                  Jindx=counteq1[(ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+ltime].setsh[0])*ha_set[arSet[ha_eqreg[eqindx]].setid].size+ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+lreg].setsh[0]]+rowindx;
                else {
                  for(i4=1; i4<MAXSUPSET; i4++)if(ha_set[arSet[ha_eqtime[eqindx]].setid].subsetid[i4]=alltimeset)break;
                  Jindx=counteq1[(ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+ltime].setsh[i4])*ha_set[arSet[ha_eqreg[eqindx]].setid].size+ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+lreg].setsh[0]]+rowindx;
                }
                //Jindx=counteq1[(ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+ltime].setsh[0])*ha_set[arSet[ha_eqreg[eqindx]].setid].size+ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+lreg].setsh[0]]+rowindx;
              }
            }
            //if(lj==0){printf("OKKKK4 istart %d iend %d j %d\n",Istart,Iend,Jindx);printf("OK\n");}
            if(alltimeset>=0&&allregset<0) {
              if(arSet[ha_eqtime[eqindx]].setid==alltimeset)
                Jindx=counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+ltime].setsh[0]]+rowindx;
              else {
                for(i4=1; i4<MAXSUPSET; i4++)if(ha_set[arSet[ha_eqtime[eqindx]].setid].subsetid[i4]=alltimeset)break;
                Jindx=counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+ltime].setsh[i4]]+rowindx;
              }
            }
            if(alltimeset<0&&allregset>=0) Jindx=counteq1[ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+ltime].setsh[0]]+rowindx;
          }
          else Jindx=matrow+rowindx;
          leq=matroworg+(uvadd)lj/nloopsfac;//matroworg+rowindxorg;
          //if(Jindx==22439)printf("line %s\n leq %d jindx %d Indx %d\n",line,leq,Jindx,Iindx);
          ha_eqadd[leq]=Jindx;
          //Jindx=matrow+rowindx;
          Iindx=ha_cgeshock[ha_var[LinVars[i].LinVarIndx].begadd+li3].ExoIndx;
          //if(Jindx==26422667)printf("name %s ltime %d rowindx %d Jindx %d Iindx %d lj %d timeindx %d fdimlin %d set %d\n",ha_var[LinVars[i].LinVarIndx].cofname,ltime,rowindx,Jindx,Iindx,lj,ha_eqtime[eqindx],fdimlin,ha_setele[arSet[ha_eqtime[eqindx]].SetBegAdd+ltime].setsh);
          //if(Jindx==26422667)for (dcount=0; dcount<fdim; dcount++)printf("dcount %d\n",dcountdim3[dcount]);
          //if(Jindx==26422667)printf("dnz %d shock %d Iend %d\n",dnnz[Jindx-Istart],ha_cgeshock[ha_var[LinVars[i].LinVarIndx].begadd+li3].ShockId,Iend);
          //Jindx=0;
          if(Istart<=Jindx&&Jindx<Iend) {
            if (!ha_cgeshock[ha_var[LinVars[i].LinVarIndx].begadd+li3].ShockId) {
              if (Istart<=Iindx&&Iindx<Iend) {
                dnnz[Jindx-Istart]=dnnz[Jindx-Istart]+1;
              }
              else {
                onnz[Jindx-Istart]=onnz[Jindx-Istart]+1;
              }
            }
            if (ha_cgeshock[ha_var[LinVars[i].LinVarIndx].begadd+li3].ShockId) {
              //PetscSynchronizedPrintf(PETSC_COMM_SELF,"matrow %d\n",Iindx);
              if (Istart<=Iindx&&Iindx<Iend) { //&&Iindx<nexo
                dnnzB[Jindx-Istart]=dnnzB[Jindx-Istart]+1;
              }
              else {
                onnzB[Jindx-Istart]=onnzB[Jindx-Istart]+1;
              }
            }
          }
          //if(lj==0)printf("OKKKK4\n");
          //if(Jindx==26422667)printf("dnz %d\n",dnnz[Jindx-Istart]);
        }
      }
      //printf("OKKK\n");
      //if(Istart<=matrow&&matrow<Iend) PetscSynchronizedPrintf(PETSC_COMM_SELF,"matrow %d dnnzA %d onnzA %d dnnzB %d onnzB %d Is %d Ie %d Jin %d\n",matrow,dnnz[matrow-Istart],onnz[matrow-Istart],dnnzB[matrow-Istart],onnzB[matrow-Istart],Istart,Iend,Jindx-Istart);
      if(ha_eqint[eqindx]) {
        //printf("beg %d dnz %d line: %s\n",counteq1[0],dnnz[counteq1[0]],line);
        if(alltimeset>=0&&allregset>=0) {
          if(nesteddbbd==1) {
            if(ha_eqreg[eqindx]>-1)for(lj=0; lj<ha_set[arSet[ha_eqtime[eqindx]].setid].size; lj++)for(l2=0; l2<ha_set[arSet[ha_eqreg[eqindx]].setid].size; l2++) {
                  counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+lj].setsh[ha_set[arSet[ha_eqtime[eqindx]].setid].intsup]*(nreg+1)+ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+l2].setsh[ha_set[arSet[ha_eqreg[eqindx]].setid].regsup]]+=nloops/ha_set[arSet[ha_eqtime[eqindx]].setid].size/ha_set[arSet[ha_eqreg[eqindx]].setid].size;
                  if(ha_set[arSet[ha_eqreg[eqindx]].setid].regsup!=0)printf("Sub reg not supported in NDBBD!!!\n");
                }
//                  counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+lj].setsh[ha_set[arSet[ha_eqtime[eqindx]].setid].intsup]*(nreg+1)+ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+l2].setsh[ha_set[arSet[ha_eqreg[eqindx]].setid].regsup]]+=nloops/ha_set[arSet[ha_eqtime[eqindx]].setid].size/ha_set[arSet[ha_eqreg[eqindx]].setid].size;
            else for(lj=0; lj<ha_set[arSet[ha_eqtime[eqindx]].setid].size; lj++)
                counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+lj].setsh[ha_set[arSet[ha_eqtime[eqindx]].setid].intsup]*(nreg+1)+nreg]+=nloops/ha_set[arSet[ha_eqtime[eqindx]].setid].size;
          }
          else {
            if(arSet[ha_eqtime[eqindx]].setid==alltimeset)
              for(lj=0; lj<ha_set[arSet[ha_eqtime[eqindx]].setid].size; lj++)for(l2=0; l2<ha_set[arSet[ha_eqreg[eqindx]].setid].size; l2++)
                  counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+lj].setsh[0]*ha_set[arSet[ha_eqreg[eqindx]].setid].size+ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+l2].setsh[0]]+=nloops/ha_set[arSet[ha_eqtime[eqindx]].setid].size/ha_set[arSet[ha_eqreg[eqindx]].setid].size;
            else {
              for(lj=0; lj<ha_set[arSet[ha_eqtime[eqindx]].setid].size; lj++)for(l2=0; l2<ha_set[arSet[ha_eqreg[eqindx]].setid].size; l2++) {
                  for(i4=1; i4<MAXSUPSET; i4++)if(ha_set[arSet[ha_eqtime[eqindx]].setid].subsetid[i4]=alltimeset)break;
                  counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+lj].setsh[i4]*ha_set[arSet[ha_eqreg[eqindx]].setid].size+ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+l2].setsh[0]]+=nloops/ha_set[arSet[ha_eqtime[eqindx]].setid].size/ha_set[arSet[ha_eqreg[eqindx]].setid].size;
                }
            }
            //counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+lj].setsh[0]*ha_set[arSet[ha_eqreg[eqindx]].setid].size+ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+l2].setsh[0]]+=nloops/ha_set[arSet[ha_eqtime[eqindx]].setid].size/ha_set[arSet[ha_eqreg[eqindx]].setid].size;
          }
        }
        //if(alltimeset>=0&&allregset<0)for(lj=0; lj<ha_set[arSet[ha_eqtime[eqindx]].setid].size; lj++)counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+lj].setsh[0]]+=nloops/ha_set[arSet[ha_eqtime[eqindx]].setid].size;
        if(alltimeset>=0&&allregset<0) {
          if(arSet[ha_eqtime[eqindx]].setid==alltimeset)
            for(lj=0; lj<ha_set[arSet[ha_eqtime[eqindx]].setid].size; lj++)counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+lj].setsh[0]]+=nloops/ha_set[arSet[ha_eqtime[eqindx]].setid].size;
          else {
            for(i4=1; i4<MAXSUPSET; i4++)if(ha_set[arSet[ha_eqtime[eqindx]].setid].subsetid[i4]=alltimeset)break;
            for(lj=0; lj<ha_set[arSet[ha_eqtime[eqindx]].setid].size; lj++)counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+lj].setsh[i4]]+=nloops/ha_set[arSet[ha_eqtime[eqindx]].setid].size;
          }
        }
        if(alltimeset<0&&allregset>=0)for(lj=0; lj<ha_set[arSet[ha_eqreg[eqindx]].setid].size; lj++)counteq1[ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+lj].setsh[0]]+=nloops/ha_set[arSet[ha_eqreg[eqindx]].setid].size;
      }
      else matrow+=nloops;
      matroworg+=nloops;
      //printf("matrow %d\n",matroworg);
      eqindx++;
      //matrow=matrow+nloops;
      free(LinVars);
      free(arSet);
      //printf("matrow %d\n",matroworg);
    }
  }
  //printf("HHHH\n");
  free(counteq1);
  fclose(filehandle);
  return 1;
}
/*
int HaReplVar(char *fname, char *commsyntax,ha_cgeset *ha_set,uvadd nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele,ha_cgeexovar *ha_cgeshock,Mat A,Mat B) {
  FILE * filehandle;
  char tline[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],leftline[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[TABREADLINE],sumsyntax[NAMESIZE],lintmp[TABREADLINE];//,*p1=NULL;
  uvdim fdim,totalsum,sumcount=1,sumcount1=0,np;
  uvadd i,rowindx,l,l1,i1=0,lvar,lvar1,lvar2,lvar3,lvar4;
  char *readitem=NULL,*p=NULL,*p1=NULL;//,*p2=NULL,*varpnts;
  ha_cgetype vval;
  PetscInt j=0,Iindx,Jindx,Istart,Iend;//,jfix=0
  PetscErrorCode ierr;
  ha_cgetype zerodivide=0;
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);
  CHKERRQ(ierr);
  uvadd dcount,dcountdim1[MAXVARDIM],dcountdim2[2*MAXVARDIM],nloops,nloopslin,fdimlin=0,i3,i4,nsumele,nsumele1,l2,sizelinvars,matrow;//,totlinvars,templinvars
  uvadd nlinvars,sumindx;
  uvdim npow,npar,nmul,nplu,ndiv,nmin,ha_calvarsize=0;

  filehandle = fopen(fname,"r");
  matrow=0;
  //printf("OK!!!!\n");
  strcpy(newtabfile1,newtabfile);
  ha_cgefrstr(newtabfile1,".","1.");
  fout = fopen(newtabfile1,"w");

  while (ha_cgertabl1(commsyntax,filehandle,line,&zerodivide)) {
    //printf("line %s\n",line);
    if (strstr(line,"(DEFAULT")==NULL) {
      ha_cgefrstr1(line, commsyntax, "");
      ha_cgefrstr1(line, "(LINEAR)", "");
      while (ha_cgefrstr(line,"  ", " "));
      while (ha_cgefrchr(line, '[', '('));
      while (ha_cgefrchr(line, ']', ')'));
      while (ha_cgefrchr(line, '{', '('));
      while (ha_cgefrchr(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=ha_cgenfind(line, "(all,");
      //printf("line %s\n",line);
      if (fdim==0) {
        readitem = strtok(line+1," ");
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      } else {
        i=ha_cgerevfind(line, "(all,");
        readitem=line+i;
        readitem = strtok(readitem,")");
        readitem = strtok(NULL,"=");
        //printf("left %s\n",readitem);
        strcpy(vname,readitem);
        //printf("left %s linecopy %s\n",linecopy1,linecopy);
        strcpy(line,linecopy);
        //printf("left %s\n",line);
        readitem = strtok(line,"=");
        //printf("read %s\n",readitem);
        readitem = strtok(NULL,";");
        //printf("left %s\n",readitem);
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      while (ha_cgefrstr(readitem," ", ""));
      while (ha_cgerecovar(readitem)==1);
      strcpy(tline,readitem);
      //printf("row %d: %s\n",matrow,readitem);
      strcpy(line1,readitem);
      npow=ha_cgenchf(readitem, '^');
      nmul=ha_cgenchf(readitem, '*');
      ndiv=ha_cgenchf(readitem, '/');
      nmul=nmul+ndiv;
      nplu=ha_cgenchf(readitem, '+');
      nmin=ha_cgenchf(readitem, '-');
      nplu=nplu+nmin;
      npar=ha_cgenchf(readitem, '(');

      strcpy(readitem,line1);
      fprintf(fout,"%s\n",readitem);    }
  }
  fclose(filehandle);
  fclose(fout);
  return 1;
}
*/

int HaDBBDMatOder(Mat A, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend, uvadd nvarele, uvadd *ha_eqadd,int *ha_rows,int *ha_cols, uvadd ndblock,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,ha_cgetype cntl6) {
  IS *rowindices,*colindices;//,isrow,iscol;
  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nrowcolin,sumrowcolin;
  Mat *submatA;
  //MPI_Status   status;
  //Mat_SeqAIJ *aa;
  PetscInt i,j,j0,j1,j2,j3,j4,nrow,ncol,nz,nz1,*ai,*aj,la;
  PetscScalar *vals;
  PetscErrorCode ierr;
  PetscViewer viewer;
  uvadd lasize;
  int *ha_rows1= (int *) calloc (VecSize,sizeof(int));
  int *ha_cols1= (int *) calloc (VecSize,sizeof(int));
  int *ha_ndblocks1= (int *) calloc (ndblock,sizeof(int));

  //*****New order test**********//
  nmatin=(uvadd)ndblock/mpisize;
  if(nmatin*mpisize<ndblock)nmatinplus=nmatin+1;
  else nmatinplus=nmatin;
  for(i=0; i<mpisize; i++)if(rank+1<=ndblock-mpisize*nmatin)nmatin++;
  //PetscInt *ha_varadds= (PetscInt *) calloc (VecSize,sizeof(PetscInt));
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  ierr = PetscMalloc(nmatin*sizeof(IS **),&rowindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&colindices);
  CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatA);
  //CHKERRQ(ierr);
  printf("rank %d nmatin %d\n",rank,nmatin);
  begblock[rank]=nmatin;
  for(i=0; i<mpisize; i++) {
    j=rank;
    //j2=j;
    MPI_Bcast(&j, 1, MPI_INT,i, PETSC_COMM_WORLD);
    //j=j2;
    MPI_Bcast(&begblock[j], 1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  j2=begblock[0];
  begblock[0]=0;
  for(i=1; i<mpisize; i++) {
    j1=begblock[i];
    begblock[i]=begblock[i-1]+j2;
    j2=j1;
  }
  //j2=begblock[0];
  //for(i=0; i<mpisize; i++)begblock[i]-=j2;
  //for(i=0; i<mpisize; i++)printf("rank %d i %d begblock %d\n",rank,i,begblock[i]);
  //ierr = PetscMalloc(nmatin*sizeof(Mat **),&submatA);CHKERRQ(ierr);
  //for(i=0;i<ndblock;i++)
  //printf("rank %d couneq %d noadd %d\n",rank,counteq[begblock[rank]],counteqnoadd[begblock[rank]]);
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    bend=counteqnoadd[i+begblock[rank]]+counteq[i+begblock[rank]];
    //if(rank==1)printf("rank %d bf %d be %d\n",rank,bfirst,bend);
//    PetscInt *indices= (PetscInt *) calloc (bend-bfirst,sizeof(PetscInt));
//    for(j=0; j<bend-bfirst; j++) {
//      indices[j]=bfirst+j;
//      //if(rank==0&&i==1)printf("indices %d rank %d\n",indices[j],ha_ndblocks[i+rank*nmatinplus]);
//    }
    //printf("rank %d bfirst %d bend %d cuneq %d\n",rank,bfirst,bend,i+begblock[rank]);
    /*if(bfirst<Istart&&i<nmatinplus-1) {
      bfirst=Istart;
      bend=Istart;
    }
    if(bfirst>=Iend) {
      bfirst=Iend;
      bend=Iend;
    }*/
    //if(rank==3)
    //printf("i %d Istart %d Iend %d first %d end %d mpisize %d ndblock %d rank %d\n",i+begblock[rank],Istart,Iend,bfirst,bend,mpisize,ndblock,rank);
    //if(i<nmatin)
    ISCreateStride(PETSC_COMM_SELF,bend-bfirst,bfirst,1,rowindices+i);
//    ISCreateGeneral(PETSC_COMM_SELF,bend-bfirst,indices,PETSC_COPY_VALUES,rowindices+i);
//    free(indices);
    //ISView(rowindices[i],PETSC_VIEWER_STDOUT_WORLD);
    bfirst=countvarintra1[i+begblock[rank]];
    bend=countvarintra1[i+1+begblock[rank]];
    //printf("rank %d bfirst1 %d bend1 %d cuneq1 %d\n",rank,bfirst,bend,i+begblock[rank]);
    //PetscInt *indices1= (PetscInt *) calloc (bend-bfirst,sizeof(PetscInt));
    //for(j=0; j<bend-bfirst; j++) {
    //indices1[j]=bfirst+j;
    //}
    //printf("i %d Istart %d Iend %d first %d end %d mpisize %d ndblock %d rank %d\n",i,Istart,Iend,bfirst,bend,mpisize,ndblock,rank);
    //if(i<nmatin)
    ISCreateStride(PETSC_COMM_SELF,bend-bfirst,bfirst,1,colindices+i);
    //ISCreateGeneral(PETSC_COMM_SELF,bend-bfirst,indices1,PETSC_COPY_VALUES,colindices+i);
    //free(indices1);
//    if(rank==1&&i==0){
//      printf("OKKKK bend-bfirst %d bfirst %d\n",bend-bfirst,bfirst);
//      ISView(colindices[i],PETSC_VIEWER_STDOUT_SELF);
//    }
    //}
  }
  printf("OK???? rank %d\n",rank);
  //if(rank==0)ISView(rowindices[1],PETSC_VIEWER_STDOUT_SELF);
  //for (i=0; i<(uvadd)ndblock/mpisize+1; i++){
  //  if(rank*((uvadd)ndblock/mpisize)+i<ndblock){
  //    ierr = MatCreateSubMatrices(A,nmatin,rowindices+rank*nmatinplus,colindices+rank*nmatinplus,MAT_INITIAL_MATRIX,&submatA);CHKERRQ(ierr);
  //  }
  //}
  //for(i=0;i<ndblock;i++){
  ierr = MatCreateSubMatrices(A,nmatin,rowindices,colindices,MAT_INITIAL_MATRIX,&submatA);
  //printf("OK0????\n");
  CHKERRQ(ierr);
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&rowindices[i]);
    CHKERRQ(ierr);
    //ierr = ISDestroy(colindices+i);CHKERRQ(ierr);
  }
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&colindices[i]);
    CHKERRQ(ierr);
  }
  PetscFree(rowindices);
  PetscFree(colindices);
  //}
  //printf("OK1????\n");
  //if(rank==0) {
  //PetscViewerASCIIOpen(PETSC_COMM_SELF, "Asub.txt", &viewer);
  //ierr = MatView(submatA[0],viewer);
  //CHKERRQ(ierr);
  //PetscViewerDestroy(&viewer);
  //}
  //MatCopy(subA,submatA[0],SAME_NONZERO_PATTERN);
  MPI_Barrier(PETSC_COMM_WORLD);
  //PetscInt *ha_eqordered= (PetscInt *) calloc (VecSize,sizeof(PetscInt));
  nrowcolin=0;
  for(j1=0; j1<nmatin; j1++) {
    Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatA[j1]->data;//*aa=subA->data;
    ai= aa->i;
    aj= aa->j;
    vals=aa->a;
    nz=aa->nz;
    nrow=submatA[j1]->rmap->n;
    ncol=submatA[j1]->cmap->n;
    printf("rank %d nz %d\n",rank,nz);
    nz1=nz;
    if(nz1<nrow)nz1=nrow;
    if(nz1<ncol)nz1=ncol;
    //la=3;
    lasize=ceil((laA/100.0)*nz1);
    int *irn=(int *) calloc (lasize,sizeof(int));
    int *jcn=(int *) calloc (lasize,sizeof(int));
    int *irn1=(int *)malloc(nrow*sizeof(int));
    int *jcn1=(int *)malloc(ncol*sizeof(int));
    ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
    int *keep=(int *) malloc ((nrow+5*ncol+4*ncol+7)*sizeof(int));//ICNTL(6)=1
  ha_cgetype *w51;
  if(nrow>ncol) w51=(ha_cgetype *) malloc ((5*nrow)*sizeof(ha_cgetype));
  else w51=(ha_cgetype *) malloc ((5*ncol)*sizeof(ha_cgetype));
  int *iw51=(int *) malloc ((6*nrow+3*ncol)*sizeof(int));
    int *insize=(int *) calloc (6,sizeof(int));
    insize[0]=nrow;
    insize[1]=ncol;
    insize[2]=nz;
    insize[4]=laA;
    insize[5]=lasize;
    for(i=0; i<nrow-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
        irn[j]=i+1;
        jcn[j]=aj[j]+1;
        values[j]=vals[j];
      }
    for(j=ai[nrow-1]; j<nz; j++) {
      irn[j]=nrow;
      jcn[j]=aj[j]+1;
      values[j]=vals[j];
    }
    //printf("irn %d jcn %d\n",irn[0],jcn[0]);
    //printf("irn %d jcn %d\n",irn[1],jcn[1]);
    //printf("irn %d jcn %d\n",irn[2],jcn[2]);
    //printf("irn %d jcn %d\n",irn[3],jcn[3]);
    j=0;
    //for(i=0; i<nz; i++) {
    //jcn[i]=aj[i]+1;
    //values[i]=vals[i];
    //}
    //memcpy(values,vals,nz);
    //for(i=0; i<nz; i++){irn[i]=ai[i]+1;jcn[i]=aj[i]+1;values[i]=vals[i];}
    //for(i=0; i<10; i++) printf("x %f irn %d jcn %d n %d m %d\n",values[i],irn[i],jcn[i],nrow,ncol);
    //printf("OK1????\n");
    spec51m_rank_(insize,&cntl6,irn,jcn,values,irn1,jcn1,keep,w51,iw51);
    //printf("OK????\n");
    for(i=0; i<nrow; i++) {
      ha_rows1[i+counteq[j1+begblock[rank]]]=irn1[i]-1;
      //printf("rank %d row %d\n",rank,ha_rows[i+counteq[j1+begblock[rank]]]);
      //if(rank==0&&j1==)printf("x %f irn %d jcn %d n %d m %d rank %d\n",values[i],irn[i]-1,jcn[i]-1,nrow,ncol,insize[3]);
    }
    for(i=0; i<ncol; i++) {
      ha_cols1[i+countvarintra1[j1+begblock[rank]]]=jcn1[i]-1;
      //if(rank==0)printf("rank %d col %d\n",rank,ha_cols[i+countvarintra1[j1+begblock[rank]]]);
      //if(rank==0)printf("rank %d i %d ha_col %d\n",rank,i,ha_cols[i+countvarintra1[j1+begblock[rank]]]);
    }
    printf("rank %d j1 %d proc %d\n",insize[3],j1,rank);
    ha_ndblocks1[j1+begblock[rank]]=insize[3];
    free(insize);
    free(irn);
    free(jcn);
    free(irn1);
    free(jcn1);
    free(values);
    free(keep);
    free(w51);
    free(iw51);
    MatDestroy(&submatA[j1]);
    nrowcolin+=ncol;
  }
  PetscFree(submatA);
  MPI_Barrier(PETSC_COMM_WORLD);
  /*j1=0;
  j=0;
  for(i=0; i<mpisize; i++) {
    j=nrowcolin;
    //j2=j;
    MPI_Bcast(&j, 1, MPI_INT,i, PETSC_COMM_WORLD);
    //j=j2;
    MPI_Bcast(&ha_cols[j1], j, MPI_INT,i, PETSC_COMM_WORLD);
    j1+=nrowcolin;
    //j2=j1;
    MPI_Bcast(&j1,1, MPI_INT,i, PETSC_COMM_WORLD);
    //j1=j2;
  }*/
  MPI_Allreduce(ha_cols1,ha_cols,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  free(ha_cols1);
  MPI_Allreduce(ha_rows1,ha_rows,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  free(ha_rows1);
  MPI_Allreduce(ha_ndblocks1,ha_ndblocks,ndblock, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  free(ha_ndblocks1);
  //for (i=0;i<ndblock;i++)printf("block %d\n",ha_ndblocks[i]);
  /*j1=0;
  j=0;
  MPI_Barrier(PETSC_COMM_WORLD);
  for(i=0; i<mpisize; i++) {
    //printf("nmatin %d proc %d rank %d j1 %d\n",nmatin,i,ha_ndblocks[j1],j1);
    j=nmatin;
    //j2=j;
    MPI_Bcast(&j,1, MPI_INT,i, PETSC_COMM_WORLD);
    //j=j2;
    MPI_Bcast(&ha_ndblocks[j1], j, MPI_INT,i, PETSC_COMM_WORLD);
    j1+=nmatin;
    //j2=j1;
    MPI_Bcast(&j1,1, MPI_INT,i, PETSC_COMM_WORLD);
    //j1=j2;
  }
  nrowcolin=0;
  for(i=0; i<nmatin; i++) {
    nrowcolin+=counteqnoadd[i+begblock[rank]];
  }
  j1=0;
  j=0;
  MPI_Barrier(PETSC_COMM_WORLD);
  printf("OK1234 rank %d\n",rank);
  for(i=0; i<mpisize; i++) {
    j=nrowcolin;
    //j2=j;
    MPI_Bcast(&j,1, MPI_INT,i, PETSC_COMM_WORLD);
    //j=j2;
    MPI_Bcast(&ha_rows[j1], j, MPI_INT,i, PETSC_COMM_WORLD);
    j1+=nrowcolin;
    //j2=j1;
    MPI_Bcast(&j1,1, MPI_INT,i, PETSC_COMM_WORLD);
    //j1=j2;
  }*/
  //free(ha_eqordered);
  //for(i=0;i<VecSize;i++)printf("rank %d col %d\n",rank,ha_cols[i]);
  free(begblock);
  MPI_Barrier(PETSC_COMM_WORLD);
  return 1;
}

//int HaDBBDMatOder1(Mat A, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend, uvadd nvarele, uvadd *ha_eqadd,int *ha_rows,int *ha_cols, uvadd ndblock,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,ha_cgetype cntl6) {
//  IS *rowindices,*colindices;//,isrow,iscol;
//  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nrowcolin,sumrowcolin;
//  Mat *submatA;
//  PetscInt i,j,j0,j1,j2,j3,j4,nrow,ncol,nz,nz1,*ai,*aj,la;
//  PetscScalar *vals;
//  PetscErrorCode ierr;
//  PetscViewer viewer;
//  int *ha_rows1= (int *) calloc (VecSize,sizeof(int));
//  int *ha_cols1= (int *) calloc (VecSize,sizeof(int));
//  int *ha_ndblocks1= (int *) calloc (ndblock,sizeof(int));
//
//  //*****New order test**********//
//  nmatin=(uvadd)ndblock/mpisize;
//  if(nmatin*mpisize<ndblock)nmatinplus=nmatin+1;
//  else nmatinplus=nmatin;
//  for(i=0; i<mpisize; i++)if(rank+1<=ndblock-mpisize*nmatin)nmatin++;
//  //PetscInt *ha_varadds= (PetscInt *) calloc (VecSize,sizeof(PetscInt));
//  int *begblock= (int *) calloc (mpisize,sizeof(int));
//  ierr = PetscMalloc(nmatin*sizeof(IS **),&rowindices);
//  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatin*sizeof(IS **),&colindices);
//  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatA);
//  CHKERRQ(ierr);
//  printf("rank %d nmatin %d\n",rank,nmatin);
//  begblock[rank]=nmatin;
//  for(i=0; i<mpisize; i++) {
//    j=rank;
//    MPI_Bcast(&j, 1, MPI_INT,i, PETSC_COMM_WORLD);
//    MPI_Bcast(&begblock[j], 1, MPI_INT,i, PETSC_COMM_WORLD);
//  }
//  j2=begblock[0];
//  begblock[0]=0;
//  for(i=1; i<mpisize; i++) {
//    j1=begblock[i];
//    begblock[i]=begblock[i-1]+j2;
//    j2=j1;
//  }
//  for(i=0; i<nmatin; i++) {
//    bfirst=counteq[i+begblock[rank]];
//    bend=counteqnoadd[i+begblock[rank]]+counteq[i+begblock[rank]];
//    ISCreateStride(PETSC_COMM_SELF,bend-bfirst,bfirst,1,rowindices+i);
//    bfirst=countvarintra1[i+begblock[rank]];
//    bend=countvarintra1[i+1+begblock[rank]];
//    ISCreateStride(PETSC_COMM_SELF,bend-bfirst,bfirst,1,colindices+i);
//  }
//  printf("OK???? rank %d\n",rank);
//  ierr = MatCreateSubMatrices(A,nmatin,rowindices,colindices,MAT_INITIAL_MATRIX,&submatA);
//  CHKERRQ(ierr);
//  for (i=0; i<nmatin; i++) {
//    ierr = ISDestroy(&rowindices[i]);
//    CHKERRQ(ierr);
//  }
//  for (i=0; i<nmatin; i++) {
//    ierr = ISDestroy(&colindices[i]);
//    CHKERRQ(ierr);
//  }
//  PetscFree(rowindices);
//  PetscFree(colindices);
//  MPI_Barrier(PETSC_COMM_WORLD);
//  nrowcolin=0;
//  for(j1=0; j1<nmatin; j1++) {
//    Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatA[j1]->data;//*aa=subA->data;
//    ai= aa->i;
//    aj= aa->j;
//    vals=aa->a;
//    nz=aa->nz;
//    nrow=submatA[j1]->rmap->n;
//    ncol=submatA[j1]->cmap->n;
//    printf("rank %d nz %d\n",rank,nz);
//    nz1=nz;
//    if(nz1<nrow)nz1=nrow;
//    if(nz1<ncol)nz1=ncol;
//    //la=3;
//    int *irn=(int *) calloc (laA*nz1,sizeof(int));
//    int *jcn=(int *) calloc (laA*nz1,sizeof(int));
//    ha_cgetype *values= (ha_cgetype *) calloc (laA*nz1,sizeof(ha_cgetype));
//    int *insize=(int *) calloc (5,sizeof(int));
//    insize[0]=nrow;
//    insize[1]=ncol;
//    insize[2]=nz;
//    insize[4]=laA;
//    for(i=0; i<nrow-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
//        irn[j]=i+1;
//        jcn[j]=aj[j]+1;
//        values[j]=vals[j];
//      }
//    for(j=ai[nrow-1]; j<nz; j++) {
//      irn[j]=nrow;
//      jcn[j]=aj[j]+1;
//      values[j]=vals[j];
//    }
//    spec51_rank_(insize,&cntl6,irn,jcn,values);
//    for(i=0; i<nrow; i++) {
//      ha_rows1[i+counteq[j1+begblock[rank]]]=irn[i]-1+counteq[j1+begblock[rank]];
//    }
//    for(i=0; i<ncol; i++) {
//      ha_cols1[i+countvarintra1[j1+begblock[rank]]]=jcn[i]-1+countvarintra1[j1+begblock[rank]];
//    }
//    printf("rank %d j1 %d proc %d\n",insize[3],j1,rank);
//    ha_ndblocks1[j1+begblock[rank]]=insize[3];
//    free(insize);
//    free(irn);
//    free(jcn);
//    free(values);
//    MatDestroy(&submatA[j1]);
//    nrowcolin+=ncol;
//  }
//  PetscFree(submatA);
//  MPI_Barrier(PETSC_COMM_WORLD);
//  MPI_Allreduce(ha_cols1,ha_cols,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
//  free(ha_cols1);
//  MPI_Allreduce(ha_rows1,ha_rows,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
//  free(ha_rows1);
//  MPI_Allreduce(ha_ndblocks1,ha_ndblocks,ndblock, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
//  free(ha_ndblocks1);
//  free(begblock);
//  MPI_Barrier(PETSC_COMM_WORLD);
//  return 1;
//}

int HaNDBBDMatOderPre(Mat A, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int nreg, int ntime, uvadd nvarele, uvadd *ha_eqadd,int *ha_rows,int *ha_cols, uvadd ndblock,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,ha_cgetype cntl6,PetscInt* ndbbdrank,PetscBool presol) {
  FILE *presolfile;
  char j1name[1024],filename[1024],rankname[1024];
  size_t frd;
  IS *rowindices=NULL,*colindices=NULL,*rowindicesD=NULL,*colindicesD=NULL;//,isrow,iscol;
  PetscInt bfirst,bend,sblockin,nmatin,nmatint,nmatinplus,nrowcolin,sumrowcolin,rfirst,cfirst;
  Mat *submatA=NULL;//,*submatD=NULL;
  PetscInt i,j,j0,j1,j2,j3,j4,j5,j6,j7,j8,nrow,ncol,nz,nz1,*ai,*aj,la;
  uvadd lasize;
  ha_cgetype cntl6in;
  PetscScalar *vals;
  PetscErrorCode ierr;
  PetscViewer viewer;
  //const PetscInt *nindices;
  MatInfo           matinfo;
  MatGetInfo(A,MAT_LOCAL,&matinfo);
  printf("rank %d matinfo.nz_used %g\n",rank,matinfo.nz_used);
  int *ha_ndblocks1= (int *) calloc (ndblock,sizeof(int));
  //uvadd *counteq1= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
  uvadd *counteq2= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
  uvadd *counteqnoadd1= (uvadd *) calloc (ndblock,sizeof(uvadd));
  //uvadd *counteq2= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
  uvadd *counteqnoadd2= (uvadd *) calloc (ndblock,sizeof(uvadd));
  uvadd *countvarintra2= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
  //uvadd *countvarintra3= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
  uvadd *countvarintra4= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
  uvadd *countvarintra6= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
  if(rank<10)strcpy(rankname,"000");
  if(rank<100&&rank>9)strcpy(rankname,"00");
  if(rank<1000&&rank>99)strcpy(rankname,"0");
  if(rank>=1000)rankname[0]='\0';
  sprintf(filename, "%d",rank);
  strcat(rankname,filename);
  //*****New order test**********//
  nmatint=(uvadd)ntime/mpisize;
  //nmatminust=nmatint;
  //if(nmatint*mpisize<ntime)nmatinplust=nmatint+1;
  //else nmatinplust=nmatint;
  for(i=0; i<mpisize; i++)if(rank+1<=ntime-mpisize*nmatint)nmatint++;
  nmatin=(nreg+1)*nmatint;
  //nmatinplus=(nreg+1)*nmatinplust;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  begblock[rank]=nmatin;
  printf("rank %d nmatin %d nmatint %d ndblock %ld\n",rank,nmatin,nmatint,ndblock);
  for(i=0; i<mpisize; i++) {
    j=rank;
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    j=j2;
    MPI_Bcast(&begblock[j2],1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  j2=begblock[0];
  begblock[0]=0;
  for(i=1; i<mpisize; i++) {
    j1=begblock[i];
    begblock[i]=begblock[i-1]+j2;
    j2=j1;
  }
  //for(i=0; i<mpisize; i++)printf("rank %d beg block %d\n",rank,begblock[i]);
//  Mat** submatBij = (Mat**)calloc(nmatin,sizeof(Mat*));
//  for (i=0; i<nmatin; i++) {
//    submatBij[i] = (Mat*)calloc(1,sizeof(Mat));
//  }
//  Mat** submatCij = (Mat**)calloc(nmatin,sizeof(Mat*));
//  for (i=0; i<nmatin; i++) {
//    submatCij[i] = (Mat*)calloc(1,sizeof(Mat));
//  }
  ierr = PetscMalloc(nmatint*nreg*sizeof(IS **),&rowindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*nreg*sizeof(IS **),&colindices);
  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*sizeof(IS **),&rowindicesD);
//  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesD);
//  CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatint*nreg*sizeof(Mat *),&submatA);
  //CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*sizeof(Mat *),&submatD);
//  CHKERRQ(ierr);
  j2=0;
  for(i=0; i<nmatint; i++) {
    for(j=0; j<nreg; j++) {
      j1=i*(nreg+1)+j;
      bfirst=counteq[j1+begblock[rank]];
      bend=counteqnoadd[j1+begblock[rank]]+counteq[j1+begblock[rank]];
      //if(rank==1)printf("rank %d j1 %d rowb %d rowe %d\n",rank,j1,bfirst,bend);
      ISCreateStride(PETSC_COMM_SELF,bend-bfirst,bfirst,1,rowindices+j2);
      bfirst=countvarintra1[j1+begblock[rank]];
      bend=countvarintra1[j1+1+begblock[rank]];
      //if(rank==1)printf("rank %d j1 %d colb %d cole %d\n",rank,j1,bfirst,bend);
      ISCreateStride(PETSC_COMM_SELF,bend-bfirst,bfirst,1,colindices+j2);
      //printf("rank %d j2 %d\n",rank,j2);
      j2++;
    }
  }
  //printf("OK???? rank %d j2 %d\n",rank,j2);
  ierr = MatCreateSubMatrices(A,nmatint*nreg,rowindices,colindices,MAT_INITIAL_MATRIX,&submatA);
  CHKERRQ(ierr);
  j2=0;
  for (i=0; i<nmatint; i++) {
    for(j=0; j<nreg; j++) {
      ierr = ISDestroy(&rowindices[j2]);
      CHKERRQ(ierr);
      ierr = ISDestroy(&colindices[j2]);
      CHKERRQ(ierr);
      j2++;
    }
  }
  PetscFree(rowindices);
  PetscFree(colindices);
  //printf("OK???? rank %d\n",rank);
  MPI_Barrier(PETSC_COMM_WORLD);
  //printf("OK1???? rank %d\n",rank);
  //nrowcolin=0;
  int *ha_rows2= (int *) calloc (VecSize,sizeof(int));
  int *ha_cols2= (int *) calloc (VecSize,sizeof(int));
  if(cntl6==0&&SORD==0)cntl6in=1e-5;
  else cntl6in=cntl6;
  #pragma omp parallel private(j1,j2,j3,j4,j5,j6,j7,j8,ai,aj,vals,nz,nrow,ncol,nz1,i,j,bfirst,lasize) shared(cntl6in,ha_rows2,ha_cols2,ha_ndblocks1,counteqnoadd1,countvarintra2,countvarintra1,submatA)
  {
  //j8=0;
  long int lasizemax=0,nrowmax=0,ncolmax=0;
  #pragma omp for schedule (static)
  for(j1=0; j1<nmatint; j1++) {
    j8=j1*nreg;
    for(j2=0; j2<nreg; j2++) {
      Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatA[j8]->data;//*aa=subA->data;
      nz=aa->nz;
      nrow=submatA[j8]->rmap->n;
      if(nrowmax<nrow)nrowmax=nrow;
      ncol=submatA[j8]->cmap->n;
      if(ncolmax<ncol)ncolmax=ncol;
      //printf("rank %d j1 %d j2 %d j8 %d nz %d\n",rank,j1,j2,j8,nz);
      nz1=nz;
      if(nz1<nrow)nz1=nrow;
      if(nz1<ncol)nz1=ncol;
      lasize=ceil((laA/100.0)*nz1);
      //printf("rank %d j1 %d j2 %d nz1 %d lasize %ld\n",rank,j1,j2,nz1,lasize);
      if(lasizemax<lasize)lasizemax=lasize;
      j8++;
    }
  }
  lasizemax+=10;
    int *irn=(int *) calloc (lasizemax,sizeof(int));
    int *jcn=(int *) calloc (lasizemax,sizeof(int));
    int*irn1=NULL,*jcn1=NULL;
  irn1=realloc(irn1,nrowmax*sizeof(int));
  jcn1=realloc(jcn1,ncolmax*sizeof(int));
    ha_cgetype *values= (ha_cgetype *) calloc (lasizemax,sizeof(ha_cgetype));
  int *keep=(int *) malloc ((nrowmax+5*ncolmax+4*ncolmax+7)*sizeof(int));//ICNTL(6)=1
  ha_cgetype *w51;
  if(nrowmax>ncolmax) w51=(ha_cgetype *) malloc ((5*nrowmax)*sizeof(ha_cgetype));
  else w51=(ha_cgetype *) malloc ((5*ncolmax)*sizeof(ha_cgetype));
  int *iw51=(int *) malloc ((6*nrowmax+3*ncolmax)*sizeof(int));
    int *insize=(int *) calloc (6,sizeof(int));
  #pragma omp for schedule (static)
  for(j1=0; j1<nmatint; j1++) {
    j4=0;
    j5=0;
    j6=0;
    j7=0;
    j8=j1*nreg;
    for(j2=0; j2<nreg; j2++) {
      j3=j1*(nreg+1)+j2;
      //printf("rank %d j1 %d j2 %d j8 %d\n",rank,j1,j2,j8);
      Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatA[j8]->data;//*aa=subA->data;
      //printf("rankxx %d j1 %d j2 %d j8 %d\n",rank,j1,j2,j8);
      ai= aa->i;
      aj= aa->j;
      vals=aa->a;
      nz=aa->nz;
      nrow=submatA[j8]->rmap->n;
      ncol=submatA[j8]->cmap->n;
      //printf("rank %d j1 %d j2 %d j8 %d nz %d\n",rank,j1,j2,j8,nz);
      nz1=nz;
      if(nz1<nrow)nz1=nrow;
      if(nz1<ncol)nz1=ncol;
      //la=3;
      lasize=ceil((laA/100.0)*nz1);
//       int *irn=(int *) calloc (lasize,sizeof(int));
//       int *jcn=(int *) calloc (lasize,sizeof(int));
//       ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
//       int *insize=(int *) calloc (5,sizeof(int));
      insize[0]=nrow;
      insize[1]=ncol;
      insize[2]=nz;
      insize[4]=laA;
      insize[5]=lasize;
      for(i=0; i<nrow-1; i++) {
        //if(rank==0&&j8==104)printf("rank %d i %d nrow %d\n",rank,i,ai[i+1]-ai[i]);
        for(j=ai[i]; j<ai[i+1]; j++) {
          irn[j]=i+1;
          jcn[j]=aj[j]+1;
          values[j]=vals[j];
          //printf("irn %d jcn %d v %lf\n",irn[j],jcn[j],vals[j]);
        }
      }
      for(j=ai[nrow-1]; j<nz; j++) {
        irn[j]=nrow;
        jcn[j]=aj[j]+1;
        values[j]=vals[j];
        //printf("irn %d jcn %d v %lf\n",irn[j],jcn[j],vals[j]);
      }
      //if(j1==0&&j2==0)return 0;
      //j=0;
      //printf("rankz %d j1 %d j2 %d nz %d\n",rank,j1,j2,nz);
      //printf("size %d values %lf nrow %d ncol %d nz %d\n",sizeof(PetscReal),values[nz-1],jcn[nz-1],irn[nz-1],nz);
      spec51m_rank_(insize,&cntl6in,irn,jcn,values,irn1,jcn1,keep,w51,iw51);
      //if(j1==0&&j2==0)return 0;
      //printf("ranka %d j1 %d j2 %d nz %d\n",rank,j1,j2,nz);
      bfirst=counteq[j3+begblock[rank]];
      for(i=0; i<nrow; i++) {
        ha_rows2[i+counteq[j3+begblock[rank]]]=irn1[i]-1+bfirst;
        //if(rank==0&&j3<=1)printf("rank %d i %d row %d\n",rank,i+counteq2[j3+begblock[rank]],ha_rows2[i+counteq[j3+begblock[rank]]]);
      }
      bfirst=countvarintra1[j3+begblock[rank]];
      for(i=0; i<ncol; i++) {
        ha_cols2[i+countvarintra1[j3+begblock[rank]]]=jcn1[i]-1+bfirst;
        //if(rank==0&&j3==4)printf("rank %d i %d row %d\n",rank,i+countvarintra1[j3+begblock[rank]],ha_cols2[i+countvarintra1[j3+begblock[rank]]]);
      }
      printf("rank %d mat rank %d nrow %d ncol %d j3 %d proc %d\n",rank,insize[3],nrow,ncol,j3,rank);
      ha_ndblocks1[j3+begblock[rank]]=insize[3];
      counteqnoadd1[j3+begblock[rank]]=insize[3];
      countvarintra2[j3+begblock[rank]]=insize[3];
      j4+=counteqnoadd1[j3+begblock[rank]];
      j5+=countvarintra2[j3+begblock[rank]];
      j6+=counteqnoadd[j3+begblock[rank]];
      j7+=(countvarintra1[j3+1+begblock[rank]]-countvarintra1[j3+begblock[rank]]);
//       free(insize);
//       free(irn);
//       free(jcn);
//       free(values);
      MatDestroy(&submatA[j8]);
      //printf("rank1 %d mat rank %d j3 %d proc %d\n",rank,insize[3],j3,rank);
      //nrowcolin+=ncol;
      //printf("j7 %d\n",j7);
      j8++;
    }
    counteqnoadd1[j1*(nreg+1)+nreg+begblock[rank]]=j6-j4+counteqnoadd[j1*(nreg+1)+nreg+begblock[rank]];
    countvarintra2[j1*(nreg+1)+nreg+begblock[rank]]=j7-j5+countvarintra1[j1*(nreg+1)+nreg+1+begblock[rank]]-countvarintra1[j1*(nreg+1)+nreg+begblock[rank]];
    //printf("j1 %d j6 %d j4 %d eq %d eq1 %d\n",j1,j6,j4,counteqnoadd[j1*(nreg+1)+nreg+begblock[rank]],counteq[j1*(nreg+1)+nreg+1+begblock[rank]]-counteq[j1*(nreg+1)+nreg+begblock[rank]]);
  }
      free(insize);
      free(irn);
      free(jcn);
      free(irn1);
      free(jcn1);
    free(keep);
    free(w51);
    free(iw51);
      free(values);
  }
  PetscFree(submatA);
  MPI_Allreduce(counteqnoadd1,counteqnoadd2,ntime*(nreg+1), MPI_LONG, MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(countvarintra2,countvarintra6,ntime*(nreg+1), MPI_LONG, MPI_SUM,PETSC_COMM_WORLD);
//  for(i=0;i<ntime*(nreg+1);i++){
//    printf("i %d c %d c1 %d\n",i,counteqnoadd2[i],countvarintra6[i]);
//  }
  //printf("OKKKK\n");
  //j4=counteq[begblock[rank]];
  //j5=countvarintra1[begblock[rank]];
  counteq2[0]=0;
  countvarintra4[0]=0;
  for(i=1; i<ndblock; i++) {
    counteq2[i]=counteq2[i-1]+counteqnoadd2[i-1];
    countvarintra4[i]=countvarintra4[i-1]+countvarintra6[i-1];
  }
  counteq2[i]=VecSize;//Not change
  countvarintra4[i]=countvarintra4[i-1]+countvarintra6[i-1];
//  for(i=0;i<ntime*(nreg+1)+1;i++){
//    printf("i %d eq %d var %d old eq %d old var %d\n",i,counteq2[i],countvarintra4[i],counteq[i],countvarintra1[i]);
//  }

//  for(i=0; i<nmatint; i++) {
//    for(j2=0; j2<nreg; j2++) {
//      j3=i*(nreg+1)+j2;
//      //counteqnoadd2[j3+begblock[rank]]=counteqnoadd1[j3+begblock[rank]];
//      counteq1[j3+begblock[rank]]=j4;
//      j4+=counteqnoadd1[j3+begblock[rank]];
//      countvarintra3[j3+begblock[rank]]=j5;
//      j5+=countvarintra2[j3+begblock[rank]];
//    }
//    j3=i*(nreg+1)+j2;
//    //counteqnoadd2[j3+begblock[rank]]=counteqnoadd1[j3+begblock[rank]];
//    counteq1[j3+begblock[rank]]=j4;
//    j4+=counteqnoadd1[j3+begblock[rank]];
//    countvarintra3[j3+begblock[rank]]=j5;
//    j5+=countvarintra2[j3+begblock[rank]];
//  }
//  MPI_Allreduce(counteq1,counteq2,ntime*(nreg+1), MPI_LONG, MPI_SUM,PETSC_COMM_WORLD);
//  //for(i=0;i<ntime*(nreg+1)+1;i++){
//  //printf("i %d c %d c1 %d c2 %d\n",i,counteq[i],counteq1[i],counteq2[i]);
//  //}
//  counteq2[ndblock]=counteq2[ndblock-1]+counteqnoadd2[ndblock-1];
//  MPI_Allreduce(countvarintra3,countvarintra4,ntime*(nreg+1), MPI_LONG, MPI_SUM,PETSC_COMM_WORLD);
//  countvarintra4[ndblock]=countvarintra4[ndblock-1]+countvarintra6[ndblock-1];//countvarintra1[ndblock];
//  for(i=0;i<ntime*(nreg+1);i++){
//    printf("rank %d i %d c %d c1 %d c2 %d\n",rank,i,counteqnoadd2[i],countvarintra4[i],countvarintra6[i]);
//  }
  //printf("OKKKK\n");
  int *ha_rows1= (int *) calloc (VecSize,sizeof(int));
  int *ha_cols1= (int *) calloc (VecSize,sizeof(int));
  for(j1=0; j1<nmatint; j1++) {
    for(j2=0; j2<nreg; j2++) {
      j3=j1*(nreg+1)+j2;
      bfirst=counteq2[j3+begblock[rank]];
      for(i=bfirst; i<counteq2[j3+begblock[rank]]+counteqnoadd2[j3+begblock[rank]]; i++) {
        ha_rows1[i]=ha_rows2[i-bfirst+counteq[j3+begblock[rank]]];
        //if(j3==1)printf("rank %d i %d ha_cols1[i] %d\n",rank,i,ha_rows1[i]);
      }
      bfirst=countvarintra4[j3+begblock[rank]];
      for(i=bfirst; i<countvarintra4[j3+1+begblock[rank]]; i++) {
        ha_cols1[i]=ha_cols2[i-bfirst+countvarintra1[j3+begblock[rank]]];
        //if(j1==10&&j2==0)printf("rank %d i %d ha_cols1[i] %d\n",rank,i,ha_cols1[i]);
      }
    }
    j3=j1*(nreg+1)+j2;
    j4=counteq2[j3+begblock[rank]];
    j5=countvarintra4[j3+begblock[rank]];
    for(j2=0; j2<nreg; j2++) {
      j6=j1*(nreg+1)+j2;
      for(i=counteqnoadd[j6+begblock[rank]]+counteq[j6+begblock[rank]]-1; i>counteqnoadd2[j6+begblock[rank]]+counteq[j6+begblock[rank]]-1; i--) {
        ha_rows1[j4]=ha_rows2[i];
        //if(j1==0)printf("i %d ha_rows1[i] %d\n",j4,ha_rows1[j4]);
        j4++;
      }
      //printf("s %d e %d\n",countvarintra1[j6+begblock[rank]]+countvarintra2[j6+begblock[rank]],countvarintra1[j6+1+begblock[rank]]);
      for(i=countvarintra1[j6+1+begblock[rank]]-1; i>countvarintra1[j6+begblock[rank]]+countvarintra6[j6+begblock[rank]]-1; i--) {
        ha_cols1[j5]=ha_cols2[i];
        //if(j1==10&&j2==0)printf("j1 %d i %d ha_cols1[i] %d\n",j1,j5,ha_cols1[j5]);
        j5++;
      }
    }
    //bfirst=counteq2[j3+begblock[rank]];
    for(i=j4; i<j4+counteqnoadd[j3+begblock[rank]]; i++) {
      ha_rows1[i]=i-j4+counteq[j3+begblock[rank]];
      //if(j1==0)printf("i1 %d ha_rows1[i] %d\n",i,ha_rows1[i]);
    }
    //bfirst=countvarintra1[j3+begblock[rank]];
    for(i=j5; i<j5+countvarintra1[j3+1+begblock[rank]]-countvarintra1[j3+begblock[rank]]; i++) {
      ha_cols1[i]=i-j5+countvarintra1[j3+begblock[rank]];
      //if(j1==10)printf("j1 %d i1 %d ha_cols1[i] %d\n",j1,i,ha_cols1[i]);
    }
  }
  free(ha_cols2);
  free(ha_rows2);
//  printf("OKKKK111\n");
//  for(i=0; i<nmatint; i++) {
//    j3=i*(nreg+1)+nreg;
//    bfirst=counteq2[j3+begblock[rank]];
//    PetscInt *indices= (PetscInt *) calloc (counteqnoadd2[j3+begblock[rank]],sizeof(PetscInt));
//    for(j=0; j<counteqnoadd2[j3+begblock[rank]]; j++) {
//      indices[j]=ha_rows1[bfirst+j];//+bfirst;
//      //printf("i %d indices[j] %d s %d\n",i,indices[j],counteqnoadd1[j3+begblock[rank]]);
//    }
//    //printf("i %d\n",i);
//    ISCreateGeneral(PETSC_COMM_SELF,counteqnoadd2[j3+begblock[rank]],indices,PETSC_COPY_VALUES,rowindicesD+i);
//    //printf("i %d\n",i);
//    indices= realloc (indices,countvarintra6[j3+begblock[rank]]*sizeof(PetscInt));
//    bfirst=countvarintra4[j3+begblock[rank]];
//    for(j=0; j<countvarintra6[j3+begblock[rank]]; j++) {
//      indices[j]=ha_cols1[bfirst+j];//+bfirst;
//      //printf("j3 %d indices[j] %d s %d\n",j3,indices[j],j);
//    }
//    ISCreateGeneral(PETSC_COMM_SELF,countvarintra6[j3+begblock[rank]],indices,PETSC_COPY_VALUES,colindicesD+i);
//    free(indices);
//  }
//  printf("OKKKK1\n");
//  ierr = MatCreateSubMatrices(A,nmatint,rowindicesD,colindicesD,MAT_INITIAL_MATRIX,&submatD);
//  CHKERRQ(ierr);
//  printf("OKKKK2 rank %d\n",rank);
//  //MPI_Barrier(PETSC_COMM_WORLD);
//  //nrowcolin=0;
//  if(cntl6==0&&SORD==0)cntl6in=1e-5;
//  else cntl6in=cntl6;
  for(j1=0; j1<nmatint; j1++) {
    j3=j1*(nreg+1)+nreg;
    nrow=counteqnoadd2[j3+begblock[rank]];
    ncol=countvarintra6[j3+begblock[rank]];
//    Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatD[j1]->data;//*aa=subA->data;
//    ai= aa->i;
//    aj= aa->j;
//    vals=aa->a;
//    nz=aa->nz;
//    //nrow=submatD[j1]->rmap->n;
//    //ncol=submatD[j1]->cmap->n;
//    //printf("rank %d j1 %d nrow %d ncol %d nz %d\n",rank,j1,nrow,ncol,nz);
//    nz1=nz;
//    if(nz1<nrow)nz1=nrow;
//    if(nz1<ncol)nz1=ncol;
//    //la=3;
//    int *irn=(int *) calloc (laDi*nz1,sizeof(int));
//    int *jcn=(int *) calloc (laDi*nz1,sizeof(int));
//    ha_cgetype *values= (ha_cgetype *) calloc (laDi*nz1,sizeof(ha_cgetype));
//    int *insize=(int *) calloc (5,sizeof(int));
//    insize[0]=nrow;
//    insize[1]=ncol;
//    insize[2]=nz;
//    insize[4]=laDi;
//    for(i=0; i<nrow-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
//        irn[j]=i+1;
//        jcn[j]=aj[j]+1;
//        values[j]=vals[j];
//      }
//    for(j=ai[nrow-1]; j<nz; j++) {
//      irn[j]=nrow;
//      jcn[j]=aj[j]+1;
//      values[j]=vals[j];
//    }
//    j=0;
//    spec51_rank_(insize,&cntl6in,irn,jcn,values);
//    ISGetIndices(rowindicesD[j1],&nindices);
//    for(i=0; i<nrow; i++) {
//      ha_rows1[i+counteq2[j3+begblock[rank]]]=nindices[irn[i]-1];
//      //printf("i %d nindices[irn[i]-1] %d counteq[j3+begblock[rank]] %d\n",i,nindices[irn[i]-1],counteq[j3+begblock[rank]]);
//    }
//    ISRestoreIndices(rowindicesD[j1],&nindices);
//    ISGetIndices(colindicesD[j1],&nindices);
//    for(i=0; i<ncol; i++) {
//      ha_cols1[i+countvarintra4[j3+begblock[rank]]]=nindices[jcn[i]-1];
//      //if(rank==0&&j1==10)printf("i %d cols %d\n",i+countvarintra4[j3+begblock[rank]],ha_cols1[i+countvarintra4[j3+begblock[rank]]]);
//    }
//    ISRestoreIndices(colindicesD[j1],&nindices);
//    printf("rank %d mat rank %d j3 %d nrow %d ncol %d\n",rank,insize[3],j3,nrow,ncol);
//    ha_ndblocks1[j3+begblock[rank]]=insize[3];
    //j4=insize[3]+15;
    //PetscSortInt(j4,ha_cols1+countvarintra4[j3+begblock[rank]]);
    //PetscSortInt(j4,ha_rows1+counteq2[j3+begblock[rank]]);
    j4=j1+begblock[rank]/(nreg+1);
    //ndbbdrank[j4]=insize[3];

    //printf("presol %d j4 %d\n",presol,j4);
//    if(presol) {
    if(nrow<ncol)ndbbdrank[j4]=nrow;
    else ndbbdrank[j4]=ncol;
    //ndbbdrank[j4]=2500;
//    } else {
//      if(nrow<ncol)ndbbdrank[j4]=nrow;
//      else ndbbdrank[j4]=ncol;
//      //PetscSortInt(insized[0],&ha_rows1[counteq2[j3+begblock[rank]]]);
//      //PetscSortInt(insized[0],&ha_cols1[countvarintra4[j3+begblock[rank]]]);
//      int *insized=(int *) calloc (5,sizeof(int));
//      if(j1<10)strcpy(j1name,"000");
//      if(j1<100&&j1>9)strcpy(j1name,"00");
//      if(j1<1000&&j1>99)strcpy(j1name,"0");
//      if(j1>=1000)j1name[0]='\0';
//      sprintf(filename, "%d",j1);
//      strcat(j1name,filename);
//      strcpy(filename,"_rank");
//      strcat(filename,rankname);
//      strcat(filename,j1name);
//      strcat(filename,".bin");
//      //printf("filename %s\n",filename);
//      if((presolfile=fopen(filename, "r"))==NULL) {
//        printf("Cannot open file. Please run presol by setting -presol 1!\n");
//      }
//      frd=fread(insized, sizeof(int), 5, presolfile);
//      fclose(presolfile);
//      if( frd== 0)printf("File read error. No presol, take min rank!\n");
//      if(ndbbdrank[j4]>insized[3]){
//        //printf("read nrow %d ncol %d in0 %d in1 %d rank %d\n",nrow,ncol,insized[0],insized[1],insized[3]);
//        int *irn1=(int *) calloc (insized[0],sizeof(int));
//        int *jcn1=(int *) calloc (insized[1],sizeof(int));
//        strcpy(filename,"_row");
//        strcat(filename,rankname);
//        strcat(filename,j1name);
//        strcat(filename,".bin");
//        if((presolfile=fopen(filename, "r"))==NULL) {
//          printf("Cannot open file.\n");
//        }
//        frd=fread(irn1, sizeof(int), insized[0], presolfile);
//        if(frd == 0)printf("File read error.");
//        fclose(presolfile);
//        strcpy(filename,"_col");
//        strcat(filename,rankname);
//        strcat(filename,j1name);
//        strcat(filename,".bin");
//        if((presolfile=fopen(filename, "r"))==NULL) {
//          printf("Cannot open file.\n");
//        }
//        frd=fread(jcn1, sizeof(int), insized[0], presolfile);
//        if( frd== 0) printf("File read error.");
//        fclose(presolfile);
//        //if(rank==2)printf("rank %d j1 %d filename %s\n",rank,j1,filename);
//        //insized[3]=insize[3];
//        //ha_ndblocks1[j3+begblock[rank]]=insized[3];
//        //if(nrow<ncol)ndbbdrank[j4]=nrow;
//        //else ndbbdrank[j4]=ncol;
//        //for(i=0;i<ncol;i++)if(rank==0&&j1==0)printf("i %d colindx %d\n",i,ha_cols1[countvarintra4[j3+begblock[rank]]+i]);
//        int *indices= (int *) calloc (insized[0],sizeof(int));
//        for(i=0; i<insized[0]; i++) {
//          indices[i]=ha_rows1[counteq2[j3+begblock[rank]]+irn1[i]-1];
//          //if(irn1[i]-1<0||irn1[i]-1>=insized[0])printf("rowOOOKKKKOOO %d j1 %d n0 %d\n",irn1[i]-1,j1,insized[0]);
//          //if(rank==0&&j1==0)printf("i %d indi %d irn %d\n",i,indices[i],irn1[i]);
//        }
////      memcpy(&ha_rows1[counteq2[j3+begblock[rank]]],indices,insized[3]*sizeof(int));
////      //memcpy(indices,&ha_rows1[counteq2[j3+begblock[rank]]],insized[3]*sizeof(int));
////      //PetscSortInt(insized[3],indices);
////      //for(i=0;i<insized[0];i++)if(rank==0&&j1==0)printf("i %d indi %d irn %d\n",i,indices[i],irn1[i]);
//        //for(i=0;i<insized[0];i++)ha_rows1[counteq2[j3+begblock[rank]]+i]=indices[i];
//        memcpy(&ha_rows1[counteq2[j3+begblock[rank]]],indices,insized[0]*sizeof(int));
//        //PetscSortInt(nrow-insized[3],&ha_rows1[counteq2[j3+begblock[rank]]+insized[3]]);
//        //memcpy(indices,&ha_rows1[counteq2[j3+begblock[rank]]+insized[3]],(nrow-insized[3])*sizeof(int));
//        //for(i=0;i<nrow-insized[3];i++)ha_rows1[counteq2[j3+begblock[rank]]+nrow-1-i]=indices[i];
////
//        for(i=0; i<insized[0]; i++){
//          indices[i]=ha_cols1[countvarintra4[j3+begblock[rank]]+jcn1[i]-1];
//          //if(jcn1[i]-1<0||jcn1[i]-1>=insized[0])if(rank==2)printf("colOOOKKKKOOO %d i %d rank %d j1 %d n0 %d\n",jcn1[i]-1,i,rank,j1,insized[0]);
//        }
//        //for(i=0; i<insized[0]; i++)ha_cols1[countvarintra4[j3+begblock[rank]]+i]=indices[i];
////      memcpy(&ha_cols1[countvarintra4[j3+begblock[rank]]],indices,insized[3]*sizeof(int));
////      //memcpy(indices,&ha_cols1[countvarintra4[j3+begblock[rank]]],insized[3]*sizeof(int));
////      //PetscSortInt(insized[3],indices);
//        memcpy(&ha_cols1[countvarintra4[j3+begblock[rank]]],indices,insized[0]*sizeof(int));
//        //PetscSortInt(ncol-insized[3],&ha_cols1[countvarintra4[j3+begblock[rank]]+insized[3]]);
//        //memcpy(indices,&ha_cols1[countvarintra4[j3+begblock[rank]]+insized[3]],(ncol-insized[3])*sizeof(int));
//        //for(i=0;i<ncol-insized[3];i++)ha_cols1[countvarintra4[j3+begblock[rank]]+ncol-1-i]=indices[i];
//        //for(i=0;i<ncol;i++)if(rank==0&&j1==0)printf("i %d colindx %d\n",i,ha_cols1[countvarintra4[j3+begblock[rank]]+i]);
//        ndbbdrank[j4]=insized[3];
//        free(indices);
//        free(irn1);
//        free(jcn1);
//      }
//      //ndbbdrank[j4]=(int)insized[3]*0.9;
//      free(insized);
//    }

    if(ndbbdrank[j4]>0) {
      //j3=j1*(nreg+1)+nreg;
      ha_ndblocks1[j3+begblock[rank]]=ndbbdrank[j4];
      //printf("OOOOKKKKK rank %d j3 %d ndblock %d tndblock %ld\n",rank,j3+begblock[rank],ha_ndblocks1[j3+begblock[rank]],ndblock);
      //for(i=0; i<nmatint; i++) {
      //ha_ndblocks1[j3+begblock[rank]]=ndbbdrank[j4];
      //int *indices= (int *) calloc (ndbbdrank[j4],sizeof(int));
      //memcpy(indices,&ha_rows1[counteq2[j3+begblock[rank]]],ndbbdrank[j4]*sizeof(int));
      //PetscSortInt(ndbbdrank[j4],indices);
      //PetscSortInt(ndbbdrank[j4],&ha_rows1[counteq2[j3+begblock[rank]]]);
      //for(j=0; j<counteqnoadd2[j3+begblock[rank]]; j++)printf("j %d ind %d\n",j,ha_rows1[j+counteq2[j3+begblock[rank]]]);
      //memcpy(&ha_rows1[counteq2[j3+begblock[rank]]],indices,ndbbdrank[j4]*sizeof(int));
      //for(j=0; j<counteqnoadd2[j3+begblock[rank]]; j++)printf("j %d ind %d\n",j,ha_rows1[j+counteq2[j3+begblock[rank]]]);
      //memcpy(indices,&ha_cols1[countvarintra4[j3+begblock[rank]]],ndbbdrank[j4]*sizeof(int));
      //for(j=0; j<ndbbdrank; j++)printf("i %d ind %d\n",j,ha_cols1[countvarintra4[j3+begblock[rank]]+j]);
      //PetscSortInt(ndbbdrank[j4],indices);
      //PetscSortInt(ndbbdrank[j4],&ha_cols1[countvarintra4[j3+begblock[rank]]]);
      //for(j=0; j<countvarintra6[j3+begblock[rank]]; j++)printf("j4a %d i %d ind %d\n",j4,j,ha_cols1[j+countvarintra4[j3+begblock[rank]]]);
      //memcpy(&ha_cols1[countvarintra4[j3+begblock[rank]]],indices,ndbbdrank[j4]*sizeof(int));
      //for(j=countvarintra6[j3+begblock[rank]]; j<countvarintra4[j3+1+begblock[rank]]; j++)printf("j4b %d j %d ind %d\n",j4,j,ha_cols1[countvarintra4[j3+begblock[rank]]+j]);
      //free(indices);
      //}
    }

//    free(insize);
//    free(irn);
//    free(jcn);
//    free(values);
//    MatDestroy(&submatD[j1]);
    //nrowcolin+=ncol;
  }
//  for (i=0; i<nmatint; i++) {
//    ierr = ISDestroy(&rowindicesD[i]);
//    CHKERRQ(ierr);
//    ierr = ISDestroy(&colindicesD[i]);
//    CHKERRQ(ierr);
//  }
//  PetscFree(submatD);
  memcpy(counteq,counteq2,(ndblock+1)*sizeof(uvadd));
  memcpy(counteqnoadd,counteqnoadd2,(ndblock)*sizeof(uvadd));
  memcpy(countvarintra1,countvarintra4,(ndblock+1)*sizeof(uvadd));
  //for (i=0; i<ndblock; i++)printf("counteq %d ceqno %d cvar %d\n",counteq[i],counteqnoadd[i],countvarintra1[i]);
//  PetscFree(rowindicesD);
//  PetscFree(colindicesD);
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Allreduce(ha_cols1,ha_cols,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(ha_rows1,ha_rows,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(ha_ndblocks1,ha_ndblocks,ndblock, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  //for(j=0; j<counteqnoadd2[8+begblock[rank]]; j++)printf("j %d ind %d\n",j,ha_rows[j+counteq2[8+begblock[rank]]]);
  //printf("OK\n");
  free(ha_cols1);
  free(ha_rows1);
  free(ha_ndblocks1);
  free(begblock);
  //free(counteq1);
  free(counteq2);
  free(counteqnoadd1);
  free(counteqnoadd2);
  free(countvarintra2);
  //free(countvarintra3);
  free(countvarintra4);
  free(countvarintra6);
  //MPI_Barrier(PETSC_COMM_WORLD);
  return 1;
}

int HaNDBBDMatOder(Mat A, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int nreg, int ntime, uvadd nvarele, uvadd *ha_eqadd,int *ha_rows,int *ha_cols, uvadd ndblock,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,ha_cgetype cntl6,PetscInt* ndbbdrank,PetscBool presol) {
  FILE *presolfile;
  char j1name[1024],filename[1024],rankname[1024];
  size_t frd;
//  IS *rowindices,*colindices,*rowindicesD,*colindicesD;//,isrow,iscol;
  PetscInt bfirst,bend,sblockin,nmatin,nmatint,nmatinplus,nrowcolin,sumrowcolin,rfirst,cfirst;
//  Mat *submatA,*submatD;
  PetscInt i,j,j0,j1,j2,j3,j4,j5,j6,j7,j8,nrow,ncol,nz,nz1,*ai,*aj,la;
//  ha_cgetype cntl6in;
//  PetscScalar *vals;
//  PetscErrorCode ierr;
//  PetscViewer viewer;
//  const PetscInt *nindices;
//  MatInfo           matinfo;
//  MatGetInfo(A,MAT_LOCAL,&matinfo);
//  printf("rank %d matinfo.nz_used %g\n",rank,matinfo.nz_used);
//  int *ha_ndblocks1= (int *) calloc (ndblock,sizeof(int));
//  //uvadd *counteq1= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
//  uvadd *counteq2= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
//  uvadd *counteqnoadd1= (uvadd *) calloc (ndblock,sizeof(uvadd));
//  //uvadd *counteq2= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
//  uvadd *counteqnoadd2= (uvadd *) calloc (ndblock,sizeof(uvadd));
//  uvadd *countvarintra2= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
//  //uvadd *countvarintra3= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
//  uvadd *countvarintra4= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
//  uvadd *countvarintra6= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
  if(rank<10)strcpy(rankname,"000");
  if(rank<100&&rank>9)strcpy(rankname,"00");
  if(rank<1000&&rank>99)strcpy(rankname,"0");
  if(rank>=1000)rankname[0]='\0';
  sprintf(filename, "%d",rank);
  strcat(rankname,filename);
  //*****New order test**********//
  nmatint=(uvadd)ntime/mpisize;
  //nmatminust=nmatint;
  //if(nmatint*mpisize<ntime)nmatinplust=nmatint+1;
  //else nmatinplust=nmatint;
  for(i=0; i<mpisize; i++)if(rank+1<=ntime-mpisize*nmatint)nmatint++;
  nmatin=(nreg+1)*nmatint;
  //nmatinplus=(nreg+1)*nmatinplust;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  begblock[rank]=nmatin;
  printf("rank %d nmatin %d nmatint %d ndblock %ld\n",rank,nmatin,nmatint,ndblock);
  for(i=0; i<mpisize; i++) {
    j=rank;
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    j=j2;
    MPI_Bcast(&begblock[j2],1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  j2=begblock[0];
  begblock[0]=0;
  for(i=1; i<mpisize; i++) {
    j1=begblock[i];
    begblock[i]=begblock[i-1]+j2;
    j2=j1;
  }
  //for(i=0; i<mpisize; i++)printf("rank %d beg block %d\n",rank,begblock[i]);
//  Mat** submatBij = (Mat**)calloc(nmatin,sizeof(Mat*));
//  for (i=0; i<nmatin; i++) {
//    submatBij[i] = (Mat*)calloc(1,sizeof(Mat));
//  }
//  Mat** submatCij = (Mat**)calloc(nmatin,sizeof(Mat*));
//  for (i=0; i<nmatin; i++) {
//    submatCij[i] = (Mat*)calloc(1,sizeof(Mat));
//  }
//  ierr = PetscMalloc(nmatint*nreg*sizeof(IS **),&rowindices);
//  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*nreg*sizeof(IS **),&colindices);
//  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*sizeof(IS **),&rowindicesD);
//  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesD);
//  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*nreg*sizeof(Mat *),&submatA);
//  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*sizeof(Mat *),&submatD);
//  CHKERRQ(ierr);
//  j2=0;
//  for(i=0; i<nmatint; i++) {
//    for(j=0; j<nreg; j++) {
//      j1=i*(nreg+1)+j;
//      bfirst=counteq[j1+begblock[rank]];
//      bend=counteqnoadd[j1+begblock[rank]]+counteq[j1+begblock[rank]];
//      //if(rank==1)printf("rank %d j1 %d rowb %d rowe %d\n",rank,j1,bfirst,bend);
//      ISCreateStride(PETSC_COMM_SELF,bend-bfirst,bfirst,1,rowindices+j2);
//      bfirst=countvarintra1[j1+begblock[rank]];
//      bend=countvarintra1[j1+1+begblock[rank]];
//      //if(rank==1)printf("rank %d j1 %d colb %d cole %d\n",rank,j1,bfirst,bend);
//      ISCreateStride(PETSC_COMM_SELF,bend-bfirst,bfirst,1,colindices+j2);
//      //printf("rank %d j2 %d\n",rank,j2);
//      j2++;
//    }
//  }
//  //printf("OK???? rank %d j2 %d\n",rank,j2);
////  ierr = MatCreateSubMatrices(A,nmatint*nreg,rowindices,colindices,MAT_INITIAL_MATRIX,&submatA);
////  CHKERRQ(ierr);
//  j2=0;
//  for (i=0; i<nmatint; i++) {
//    for(j=0; j<nreg; j++) {
//    ierr = ISDestroy(&rowindices[j2]);
//    CHKERRQ(ierr);
//    ierr = ISDestroy(&colindices[j2]);
//    CHKERRQ(ierr);
//    j2++;
//    }
//  }
//  PetscFree(rowindices);
//  PetscFree(colindices);
//  //printf("OK???? rank %d\n",rank);
//  MPI_Barrier(PETSC_COMM_WORLD);
//  printf("OK1???? rank %d\n",rank);
//  //nrowcolin=0;
//  j8=0;
//  int *ha_rows2= (int *) calloc (VecSize,sizeof(int));
//  int *ha_cols2= (int *) calloc (VecSize,sizeof(int));
//  if(cntl6==0&&SORD==0)cntl6in=1e-5;
//  else cntl6in=cntl6;
//  for(j1=0; j1<nmatint; j1++) {
//    j4=0;
//    j5=0;
//    j6=0;
//    j7=0;
//    for(j2=0; j2<nreg; j2++) {
////      j3=j1*(nreg+1)+j2;
////      //printf("rank %d j1 %d j2 %d j8 %d\n",rank,j1,j2,j8);
////      Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatA[j8]->data;//*aa=subA->data;
////      //printf("rankxx %d j1 %d j2 %d j8 %d\n",rank,j1,j2,j8);
////      ai= aa->i;
////      aj= aa->j;
////      vals=aa->a;
////      nz=aa->nz;
////      nrow=submatA[j8]->rmap->n;
////      ncol=submatA[j8]->cmap->n;
////      printf("rank %d j1 %d j2 %d j8 %d nz %d\n",rank,j1,j2,j8,nz);
////      nz1=nz;
////      if(nz1<nrow)nz1=nrow;
////      if(nz1<ncol)nz1=ncol;
//      //la=3;
////      int *irn=(int *) calloc (laA*nz1,sizeof(int));
////      int *jcn=(int *) calloc (laA*nz1,sizeof(int));
////      ha_cgetype *values= (ha_cgetype *) calloc (laA*nz1,sizeof(ha_cgetype));
////      int *insize=(int *) calloc (5,sizeof(int));
////      insize[0]=nrow;
////      insize[1]=ncol;
////      insize[2]=nz;
////      insize[4]=laA;
////      for(i=0; i<nrow-1; i++) {
////        //if(rank==0&&j8==104)printf("rank %d i %d nrow %d\n",rank,i,ai[i+1]-ai[i]);
////        for(j=ai[i]; j<ai[i+1]; j++) {
////          irn[j]=i+1;
////          jcn[j]=aj[j]+1;
////          values[j]=vals[j];
////          //printf("irn %d jcn %d v %lf\n",irn[j],jcn[j],vals[j]);
////        }
////      }
////      for(j=ai[nrow-1]; j<nz; j++) {
////        irn[j]=nrow;
////        jcn[j]=aj[j]+1;
////        values[j]=vals[j];
////        //printf("irn %d jcn %d v %lf\n",irn[j],jcn[j],vals[j]);
////      }
////      //if(j1==0&&j2==0)return 0;
////      //j=0;
////      //printf("rankz %d j1 %d j2 %d nz %d\n",rank,j1,j2,nz);
////      //printf("size %d values %lf nrow %d ncol %d nz %d\n",sizeof(PetscReal),values[nz-1],jcn[nz-1],irn[nz-1],nz);
////      spec51_rank_(insize,&cntl6in,irn,jcn,values);
////      //if(j1==0&&j2==0)return 0;
////      //printf("ranka %d j1 %d j2 %d nz %d\n",rank,j1,j2,nz);
////      bfirst=counteq[j3+begblock[rank]];
////      for(i=0; i<nrow; i++) {
////        ha_rows2[i+counteq[j3+begblock[rank]]]=irn[i]-1+bfirst;
////        //if(rank==0&&j3<=1)printf("rank %d i %d row %d\n",rank,i+counteq2[j3+begblock[rank]],ha_rows2[i+counteq[j3+begblock[rank]]]);
////      }
////      bfirst=countvarintra1[j3+begblock[rank]];
////      for(i=0; i<ncol; i++) {
////        ha_cols2[i+countvarintra1[j3+begblock[rank]]]=jcn[i]-1+bfirst;
////        //if(rank==0&&j3==4)printf("rank %d i %d row %d\n",rank,i+countvarintra1[j3+begblock[rank]],ha_cols2[i+countvarintra1[j3+begblock[rank]]]);
////      }
////      printf("rank %d mat rank %d nrow %d ncol %d j3 %d proc %d\n",rank,insize[3],nrow,ncol,j3,rank);
////      ha_ndblocks1[j3+begblock[rank]]=insize[3];
////      counteqnoadd1[j3+begblock[rank]]=insize[3];
////      countvarintra2[j3+begblock[rank]]=insize[3];
////      j4+=counteqnoadd1[j3+begblock[rank]];
////      j5+=countvarintra2[j3+begblock[rank]];
////      j6+=counteqnoadd[j3+begblock[rank]];
////      j7+=(countvarintra1[j3+1+begblock[rank]]-countvarintra1[j3+begblock[rank]]);
////      free(insize);
////      free(irn);
////      free(jcn);
////      free(values);
////      MatDestroy(&submatA[j8]);
//      //printf("rank1 %d mat rank %d j3 %d proc %d\n",rank,insize[3],j3,rank);
//      //nrowcolin+=ncol;
//      //printf("j7 %d\n",j7);
//      j8++;
//    }
////    counteqnoadd1[j1*(nreg+1)+nreg+begblock[rank]]=j6-j4+counteqnoadd[j1*(nreg+1)+nreg+begblock[rank]];
////    countvarintra2[j1*(nreg+1)+nreg+begblock[rank]]=j7-j5+countvarintra1[j1*(nreg+1)+nreg+1+begblock[rank]]-countvarintra1[j1*(nreg+1)+nreg+begblock[rank]];
//    //printf("j1 %d j6 %d j4 %d eq %d eq1 %d\n",j1,j6,j4,counteqnoadd[j1*(nreg+1)+nreg+begblock[rank]],counteq[j1*(nreg+1)+nreg+1+begblock[rank]]-counteq[j1*(nreg+1)+nreg+begblock[rank]]);
//  }
//  PetscFree(submatA);
//  MPI_Allreduce(counteqnoadd1,counteqnoadd2,ntime*(nreg+1), MPI_LONG, MPI_SUM,PETSC_COMM_WORLD);
//  MPI_Allreduce(countvarintra2,countvarintra6,ntime*(nreg+1), MPI_LONG, MPI_SUM,PETSC_COMM_WORLD);
//  for(i=0;i<ntime*(nreg+1);i++){
//    printf("i %d c %d c1 %d\n",i,counteqnoadd2[i],countvarintra6[i]);
//  }
//  printf("OKKKK\n");
//  //j4=counteq[begblock[rank]];
//  //j5=countvarintra1[begblock[rank]];
//  counteq2[0]=0;
//  countvarintra4[0]=0;
//  for(i=1; i<ndblock; i++) {
//    counteq2[i]=counteq2[i-1]+counteqnoadd2[i-1];
//    countvarintra4[i]=countvarintra4[i-1]+countvarintra6[i-1];
//  }
//  counteq2[i]=VecSize;//Not change
//  countvarintra4[i]=countvarintra4[i-1]+countvarintra6[i-1];
//  for(i=0;i<ntime*(nreg+1)+1;i++){
//    printf("i %d eq %d var %d old eq %d old var %d\n",i,counteq2[i],countvarintra4[i],counteq[i],countvarintra1[i]);
//  }

//  for(i=0; i<nmatint; i++) {
//    for(j2=0; j2<nreg; j2++) {
//      j3=i*(nreg+1)+j2;
//      //counteqnoadd2[j3+begblock[rank]]=counteqnoadd1[j3+begblock[rank]];
//      counteq1[j3+begblock[rank]]=j4;
//      j4+=counteqnoadd1[j3+begblock[rank]];
//      countvarintra3[j3+begblock[rank]]=j5;
//      j5+=countvarintra2[j3+begblock[rank]];
//    }
//    j3=i*(nreg+1)+j2;
//    //counteqnoadd2[j3+begblock[rank]]=counteqnoadd1[j3+begblock[rank]];
//    counteq1[j3+begblock[rank]]=j4;
//    j4+=counteqnoadd1[j3+begblock[rank]];
//    countvarintra3[j3+begblock[rank]]=j5;
//    j5+=countvarintra2[j3+begblock[rank]];
//  }
//  MPI_Allreduce(counteq1,counteq2,ntime*(nreg+1), MPI_LONG, MPI_SUM,PETSC_COMM_WORLD);
//  //for(i=0;i<ntime*(nreg+1)+1;i++){
//  //printf("i %d c %d c1 %d c2 %d\n",i,counteq[i],counteq1[i],counteq2[i]);
//  //}
//  counteq2[ndblock]=counteq2[ndblock-1]+counteqnoadd2[ndblock-1];
//  MPI_Allreduce(countvarintra3,countvarintra4,ntime*(nreg+1), MPI_LONG, MPI_SUM,PETSC_COMM_WORLD);
//  countvarintra4[ndblock]=countvarintra4[ndblock-1]+countvarintra6[ndblock-1];//countvarintra1[ndblock];
//  for(i=0;i<ntime*(nreg+1);i++){
//    printf("rank %d i %d c %d c1 %d c2 %d\n",rank,i,counteqnoadd2[i],countvarintra4[i],countvarintra6[i]);
//  }
//  printf("OKKKK\n");
//  int *ha_rows1= (int *) calloc (VecSize,sizeof(int));
//  int *ha_cols1= (int *) calloc (VecSize,sizeof(int));
//  for(j1=0; j1<nmatint; j1++) {
//    for(j2=0; j2<nreg; j2++) {
//      j3=j1*(nreg+1)+j2;
//      bfirst=counteq2[j3+begblock[rank]];
//      for(i=bfirst; i<counteq2[j3+begblock[rank]]+counteqnoadd2[j3+begblock[rank]]; i++) {
//        ha_rows1[i]=ha_rows2[i-bfirst+counteq[j3+begblock[rank]]];
//        //if(j3==1)printf("rank %d i %d ha_cols1[i] %d\n",rank,i,ha_rows1[i]);
//      }
//      bfirst=countvarintra4[j3+begblock[rank]];
//      for(i=bfirst; i<countvarintra4[j3+1+begblock[rank]]; i++) {
//        ha_cols1[i]=ha_cols2[i-bfirst+countvarintra1[j3+begblock[rank]]];
//        //if(j1==10&&j2==0)printf("rank %d i %d ha_cols1[i] %d\n",rank,i,ha_cols1[i]);
//      }
//    }
//    j3=j1*(nreg+1)+j2;
//    j4=counteq2[j3+begblock[rank]];
//    j5=countvarintra4[j3+begblock[rank]];
//    for(j2=0; j2<nreg; j2++) {
//      j6=j1*(nreg+1)+j2;
//      for(i=counteqnoadd[j6+begblock[rank]]+counteq[j6+begblock[rank]]-1; i>counteqnoadd2[j6+begblock[rank]]+counteq[j6+begblock[rank]]-1; i--) {
//        ha_rows1[j4]=ha_rows2[i];
//        //if(j1==0)printf("i %d ha_rows1[i] %d\n",j4,ha_rows1[j4]);
//        j4++;
//      }
//      //printf("s %d e %d\n",countvarintra1[j6+begblock[rank]]+countvarintra2[j6+begblock[rank]],countvarintra1[j6+1+begblock[rank]]);
//      for(i=countvarintra1[j6+1+begblock[rank]]-1; i>countvarintra1[j6+begblock[rank]]+countvarintra6[j6+begblock[rank]]-1; i--) {
//        ha_cols1[j5]=ha_cols2[i];
//        //if(j1==10&&j2==0)printf("j1 %d i %d ha_cols1[i] %d\n",j1,j5,ha_cols1[j5]);
//        j5++;
//      }
//    }
//    //bfirst=counteq2[j3+begblock[rank]];
//    for(i=j4; i<j4+counteqnoadd[j3+begblock[rank]]; i++) {
//      ha_rows1[i]=i-j4+counteq[j3+begblock[rank]];
//      //if(j1==0)printf("i1 %d ha_rows1[i] %d\n",i,ha_rows1[i]);
//    }
//    //bfirst=countvarintra1[j3+begblock[rank]];
//    for(i=j5; i<j5+countvarintra1[j3+1+begblock[rank]]-countvarintra1[j3+begblock[rank]]; i++) {
//      ha_cols1[i]=i-j5+countvarintra1[j3+begblock[rank]];
//      //if(j1==10)printf("j1 %d i1 %d ha_cols1[i] %d\n",j1,i,ha_cols1[i]);
//    }
//  }
//  free(ha_cols2);
//  free(ha_rows2);
//  printf("OKKKK111\n");
//  for(i=0; i<nmatint; i++) {
//    j3=i*(nreg+1)+nreg;
//    bfirst=counteq2[j3+begblock[rank]];
//    PetscInt *indices= (PetscInt *) calloc (counteqnoadd2[j3+begblock[rank]],sizeof(PetscInt));
//    for(j=0; j<counteqnoadd2[j3+begblock[rank]]; j++) {
//      indices[j]=ha_rows1[bfirst+j];//+bfirst;
//      //printf("i %d indices[j] %d s %d\n",i,indices[j],counteqnoadd1[j3+begblock[rank]]);
//    }
//    //printf("i %d\n",i);
//    ISCreateGeneral(PETSC_COMM_SELF,counteqnoadd2[j3+begblock[rank]],indices,PETSC_COPY_VALUES,rowindicesD+i);
//    //printf("i %d\n",i);
//    indices= realloc (indices,countvarintra6[j3+begblock[rank]]*sizeof(PetscInt));
//    bfirst=countvarintra4[j3+begblock[rank]];
//    for(j=0; j<countvarintra6[j3+begblock[rank]]; j++) {
//      indices[j]=ha_cols1[bfirst+j];//+bfirst;
//      //printf("j3 %d indices[j] %d s %d\n",j3,indices[j],j);
//    }
//    ISCreateGeneral(PETSC_COMM_SELF,countvarintra6[j3+begblock[rank]],indices,PETSC_COPY_VALUES,colindicesD+i);
//    free(indices);
//  }
//  printf("OKKKK1\n");
//  ierr = MatCreateSubMatrices(A,nmatint,rowindicesD,colindicesD,MAT_INITIAL_MATRIX,&submatD);
//  CHKERRQ(ierr);
//  printf("OKKKK2 rank %d\n",rank);
//  //MPI_Barrier(PETSC_COMM_WORLD);
//  //nrowcolin=0;
//  if(cntl6==0&&SORD==0)cntl6in=1e-5;
//  else cntl6in=cntl6;
  int *ha_rows3= (int *) calloc (VecSize,sizeof(int));
  int *ha_cols3= (int *) calloc (VecSize,sizeof(int));
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    bend=counteqnoadd[i+begblock[rank]]+counteq[i+begblock[rank]];
    //printf("bfisrt %d bend %d Vs %ld\n",bfirst,bend,VecSize);
    for(j=bfirst; j<bend; j++)ha_rows3[j]=ha_rows[j];
    bfirst=countvarintra1[i+begblock[rank]];
    bend=countvarintra1[i+1+begblock[rank]];
    for(j=bfirst; j<bend; j++)ha_cols3[j]=ha_cols[j];
  }
  int *ha_ndblocks2= (int *) calloc (ndblock,sizeof(int));
  for(j1=0; j1<nmatin; j1++)ha_ndblocks2[j1+begblock[rank]]=ha_ndblocks[j1+begblock[rank]];

  for(j1=0; j1<nmatint; j1++) {
    j3=j1*(nreg+1)+nreg;
//    Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatD[j1]->data;//*aa=subA->data;
//    ai= aa->i;
//    aj= aa->j;
//    vals=aa->a;
//    nz=aa->nz;
//    nrow=submatD[j1]->rmap->n;
//    ncol=submatD[j1]->cmap->n;
//    //printf("rank %d j1 %d nrow %d ncol %d nz %d\n",rank,j1,nrow,ncol,nz);
//    nz1=nz;
//    if(nz1<nrow)nz1=nrow;
//    if(nz1<ncol)nz1=ncol;
//    //la=3;
//    int *irn=(int *) calloc (laDi*nz1,sizeof(int));
//    int *jcn=(int *) calloc (laDi*nz1,sizeof(int));
//    ha_cgetype *values= (ha_cgetype *) calloc (laDi*nz1,sizeof(ha_cgetype));
//    int *insize=(int *) calloc (5,sizeof(int));
//    insize[0]=nrow;
//    insize[1]=ncol;
//    insize[2]=nz;
//    insize[4]=laDi;
//    for(i=0; i<nrow-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
//        irn[j]=i+1;
//        jcn[j]=aj[j]+1;
//        values[j]=vals[j];
//      }
//    for(j=ai[nrow-1]; j<nz; j++) {
//      irn[j]=nrow;
//      jcn[j]=aj[j]+1;
//      values[j]=vals[j];
//    }
//    j=0;
//    //spec51_rank_(insize,&cntl6in,irn,jcn,values);
//    ISGetIndices(rowindicesD[j1],&nindices);
//    for(i=0; i<nrow; i++) {
//      //printf("bf %d aft %d\n",ha_rows3[i+counteq2[j3+begblock[rank]]],nindices[irn[i]-1]);
//      ha_rows1[i+counteq2[j3+begblock[rank]]]=nindices[irn[i]-1];
//      //printf("i %d nindices[irn[i]-1] %d counteq[j3+begblock[rank]] %d\n",i,nindices[irn[i]-1],counteq[j3+begblock[rank]]);
//    }
//    ISRestoreIndices(rowindicesD[j1],&nindices);
//    ISGetIndices(colindicesD[j1],&nindices);
//    for(i=0; i<ncol; i++) {
//      ha_cols1[i+countvarintra4[j3+begblock[rank]]]=nindices[jcn[i]-1];
//      //if(rank==0&&j1==10)printf("i %d cols %d\n",i+countvarintra4[j3+begblock[rank]],ha_cols1[i+countvarintra4[j3+begblock[rank]]]);
//    }
//    ISRestoreIndices(colindicesD[j1],&nindices);
//    printf("rank %d mat rank %d j3 %d nrow %d ncol %d\n",rank,insize[3],j3,nrow,ncol);
//    ha_ndblocks1[j3+begblock[rank]]=insize[3];
    //j4=insize[3]+15;
    //PetscSortInt(j4,ha_cols1+countvarintra4[j3+begblock[rank]]);
    //PetscSortInt(j4,ha_rows1+counteq2[j3+begblock[rank]]);
    j4=j1+begblock[rank]/(nreg+1);
    //ndbbdrank[j4]=insize[3];
    //memcpy(ha_rows1,ha_rows,VecSize*sizeof(int));
//    printf("presol %d j4 %d\n",presol,j4);
//
//      if(nrow<ncol)ndbbdrank[j4]=nrow;
//      else ndbbdrank[j4]=ncol;



    if(!presol) {
      //PetscSortInt(insized[0],&ha_rows1[counteq2[j3+begblock[rank]]]);
      //PetscSortInt(insized[0],&ha_cols1[countvarintra4[j3+begblock[rank]]]);
      int *insized=(int *) calloc (5,sizeof(int));
      if(j1<10)strcpy(j1name,"000");
      if(j1<100&&j1>9)strcpy(j1name,"00");
      if(j1<1000&&j1>99)strcpy(j1name,"0");
      if(j1>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j1);
      strcat(j1name,filename);
      strcpy(filename,"_rank");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      //printf("filename %s\n",filename);
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Cannot open file. Please run presol by setting -presol 1!\n");
      }
      frd=fread(insized, sizeof(int), 5, presolfile);
      fclose(presolfile);
      if( frd== 0)printf("File read error. No presol, take min rank!\n");
      if(ndbbdrank[j4]>insized[3]) {
        //printf("read nrow %d ncol %d in0 %d in1 %d rank %d\n",nrow,ncol,insized[0],insized[1],insized[3]);
        int *irn1=(int *) calloc (insized[0],sizeof(int));
        int *jcn1=(int *) calloc (insized[1],sizeof(int));
        strcpy(filename,"_row");
        strcat(filename,rankname);
        strcat(filename,j1name);
        strcat(filename,".bin");
        if((presolfile=fopen(filename, "r"))==NULL) {
          printf("Cannot open file.\n");
        }
        frd=fread(irn1, sizeof(int), insized[0], presolfile);
        if(frd == 0)printf("File read error.");
        fclose(presolfile);
        strcpy(filename,"_col");
        strcat(filename,rankname);
        strcat(filename,j1name);
        strcat(filename,".bin");
        if((presolfile=fopen(filename, "r"))==NULL) {
          printf("Cannot open file.\n");
        }
        frd=fread(jcn1, sizeof(int), insized[0], presolfile);
        if( frd== 0) printf("File read error.");
        fclose(presolfile);
        //if(rank==2)printf("rank %d j1 %d filename %s\n",rank,j1,filename);
        //insized[3]=insize[3];
        //ha_ndblocks1[j3+begblock[rank]]=insized[3];
        //if(nrow<ncol)ndbbdrank[j4]=nrow;
        //else ndbbdrank[j4]=ncol;
        //for(i=0;i<ncol;i++)if(rank==0&&j1==0)printf("i %d colindx %d\n",i,ha_cols1[countvarintra4[j3+begblock[rank]]+i]);
        int *indices= (int *) calloc (insized[0],sizeof(int));
        for(i=0; i<insized[0]; i++) {
          indices[i]=ha_rows3[counteq[j3+begblock[rank]]+irn1[i]-1];
          //if(irn1[i]-1<0||irn1[i]-1>=insized[0])printf("rowOOOKKKKOOO %d j1 %d n0 %d\n",irn1[i]-1,j1,insized[0]);
          //if(rank==0&&j1==0)printf("i %d indi %d irn %d\n",i,indices[i],irn1[i]);
        }
//      memcpy(&ha_rows1[counteq2[j3+begblock[rank]]],indices,insized[3]*sizeof(int));
//      //memcpy(indices,&ha_rows1[counteq2[j3+begblock[rank]]],insized[3]*sizeof(int));
//      //PetscSortInt(insized[3],indices);
//      //for(i=0;i<insized[0];i++)if(rank==0&&j1==0)printf("i %d indi %d irn %d\n",i,indices[i],irn1[i]);
        //for(i=0;i<insized[0];i++)ha_rows1[counteq2[j3+begblock[rank]]+i]=indices[i];
        memcpy(&ha_rows3[counteq[j3+begblock[rank]]],indices,insized[0]*sizeof(int));
        //PetscSortInt(nrow-insized[3],&ha_rows1[counteq2[j3+begblock[rank]]+insized[3]]);
        //memcpy(indices,&ha_rows1[counteq2[j3+begblock[rank]]+insized[3]],(nrow-insized[3])*sizeof(int));
        //for(i=0;i<nrow-insized[3];i++)ha_rows1[counteq2[j3+begblock[rank]]+nrow-1-i]=indices[i];
//
        for(i=0; i<insized[0]; i++) {
          indices[i]=ha_cols3[countvarintra1[j3+begblock[rank]]+jcn1[i]-1];
          //if(jcn1[i]-1<0||jcn1[i]-1>=insized[0])if(rank==2)printf("colOOOKKKKOOO %d i %d rank %d j1 %d n0 %d\n",jcn1[i]-1,i,rank,j1,insized[0]);
        }
        //for(i=0; i<insized[0]; i++)ha_cols1[countvarintra4[j3+begblock[rank]]+i]=indices[i];
//      memcpy(&ha_cols1[countvarintra4[j3+begblock[rank]]],indices,insized[3]*sizeof(int));
//      //memcpy(indices,&ha_cols1[countvarintra4[j3+begblock[rank]]],insized[3]*sizeof(int));
//      //PetscSortInt(insized[3],indices);
        memcpy(&ha_cols3[countvarintra1[j3+begblock[rank]]],indices,insized[0]*sizeof(int));
        //PetscSortInt(ncol-insized[3],&ha_cols1[countvarintra4[j3+begblock[rank]]+insized[3]]);
        //memcpy(indices,&ha_cols1[countvarintra4[j3+begblock[rank]]+insized[3]],(ncol-insized[3])*sizeof(int));
        //for(i=0;i<ncol-insized[3];i++)ha_cols1[countvarintra4[j3+begblock[rank]]+ncol-1-i]=indices[i];
        //for(i=0;i<ncol;i++)if(rank==0&&j1==0)printf("i %d colindx %d\n",i,ha_cols1[countvarintra4[j3+begblock[rank]]+i]);
        ndbbdrank[j4]=insized[3];
        free(indices);
        free(irn1);
        free(jcn1);
      }
      //ndbbdrank[j4]=(int)insized[3]*0.9;
      free(insized);
    }

    if(ndbbdrank[j4]>0) {
      //j3=j1*(nreg+1)+nreg;
      ha_ndblocks2[j3+begblock[rank]]=ndbbdrank[j4];
      //printf("OOOOKKKKK rank %d j3 %d ndblock %d tndblock %ld\n",rank,j3+begblock[rank],ha_ndblocks2[j3+begblock[rank]],ndblock);
      //for(i=0; i<nmatint; i++) {
      //ha_ndblocks1[j3+begblock[rank]]=ndbbdrank[j4];
      //int *indices= (int *) calloc (ndbbdrank[j4],sizeof(int));
      //memcpy(indices,&ha_rows1[counteq2[j3+begblock[rank]]],ndbbdrank[j4]*sizeof(int));
      //PetscSortInt(ndbbdrank[j4],indices);
      //PetscSortInt(ndbbdrank[j4],&ha_rows1[counteq2[j3+begblock[rank]]]);
      //for(j=0; j<counteqnoadd2[j3+begblock[rank]]; j++)printf("j %d ind %d\n",j,ha_rows1[j+counteq2[j3+begblock[rank]]]);
      //memcpy(&ha_rows1[counteq2[j3+begblock[rank]]],indices,ndbbdrank[j4]*sizeof(int));
      //for(j=0; j<counteqnoadd2[j3+begblock[rank]]; j++)printf("j %d ind %d\n",j,ha_rows1[j+counteq2[j3+begblock[rank]]]);
      //memcpy(indices,&ha_cols1[countvarintra4[j3+begblock[rank]]],ndbbdrank[j4]*sizeof(int));
      //for(j=0; j<ndbbdrank; j++)printf("i %d ind %d\n",j,ha_cols1[countvarintra4[j3+begblock[rank]]+j]);
      //PetscSortInt(ndbbdrank[j4],indices);
      //PetscSortInt(ndbbdrank[j4],&ha_cols1[countvarintra4[j3+begblock[rank]]]);
      //for(j=0; j<countvarintra6[j3+begblock[rank]]; j++)printf("j4a %d i %d ind %d\n",j4,j,ha_cols1[j+countvarintra4[j3+begblock[rank]]]);
      //memcpy(&ha_cols1[countvarintra4[j3+begblock[rank]]],indices,ndbbdrank[j4]*sizeof(int));
      //for(j=countvarintra6[j3+begblock[rank]]; j<countvarintra4[j3+1+begblock[rank]]; j++)printf("j4b %d j %d ind %d\n",j4,j,ha_cols1[countvarintra4[j3+begblock[rank]]+j]);
      //free(indices);
      //}
    }

//    free(insize);
//    free(irn);
//    free(jcn);
//    free(values);
//    MatDestroy(&submatD[j1]);
    //nrowcolin+=ncol;
  }
//  for (i=0; i<nmatint; i++) {
//    ierr = ISDestroy(&rowindicesD[i]);
//    CHKERRQ(ierr);
//    ierr = ISDestroy(&colindicesD[i]);
//    CHKERRQ(ierr);
//  }
//  PetscFree(submatD);
  //for(i=0;i<ndblock+1;i++)printf("coun %d count2 %d\n",counteq[i],counteq2[i]);
  //memcpy(counteq,counteq2,(ndblock+1)*sizeof(uvadd));
  //memcpy(counteqnoadd,counteqnoadd2,(ndblock)*sizeof(uvadd));
  //memcpy(countvarintra1,countvarintra4,(ndblock+1)*sizeof(uvadd));
  //for (i=0; i<ndblock; i++)printf("counteq %d ceqno %d cvar %d\n",counteq[i],counteqnoadd[i],countvarintra1[i]);
//  PetscFree(rowindicesD);
//  PetscFree(colindicesD);
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Allreduce(ha_cols3,ha_cols,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(ha_rows3,ha_rows,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(ha_ndblocks2,ha_ndblocks,ndblock, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  //for(j=0; j<counteqnoadd2[8+begblock[rank]]; j++)printf("j %d ind %d\n",j,ha_rows[j+counteq2[8+begblock[rank]]]);
  printf("OK\n");
//  free(ha_cols1);
//  free(ha_rows1);
  free(ha_cols3);
  free(ha_rows3);
//  free(ha_ndblocks1);
  free(ha_ndblocks2);
  free(begblock);
  //free(counteq1);
//  free(counteq2);
//  free(counteqnoadd1);
//  free(counteqnoadd2);
//  free(countvarintra2);
  //free(countvarintra3);
//  free(countvarintra4);
//  free(countvarintra6);
  //MPI_Barrier(PETSC_COMM_WORLD);
  return 1;
}

int HaDBBDSol(Mat A, Vec b, ha_cgetype *x1, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, uvadd ndblock,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laD) {
  IS *rowindices,*colindices,*Cindices,*Bindices;
  const PetscInt *nindices;
  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nrowcolin,sumrowcolin,i,j,j0,j1,j2,j3,j4,j5,j6,l0,l1,l2,l3,l4,l5,rank1,proc1=0,nnzmax,j1nz,j1irnbs;
  Mat *submatA,*submatC,*submatB,*submatD,submatBT;
  Vec vecxd;
  PetscInt *ai,*aj,*aic,*ajc;
  PetscInt nrow,ncol,nz,nrowc,ncolc,nrowb,ncolb,nzc,nzv;
  forint iv,jv,la1;
  uvadd lasize,ldsize;
  MPI_Status   status;
  clock_t timestr,timeend;
  forint forintwrp;
  int systemresult;
  size_t freadresult;
  //long long int nrow,ncol,nz,nrowc,ncolc,nzc,nzv;
  PetscScalar *vals,*valsc,vecval;
  ha_cgetype vval;
  PetscErrorCode ierr;
  PetscViewer viewer;
  FILE* fp1,*fp2,*fp3;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  char filename[1024],rankname[1024],j1name[1024],tempchar[1024];//,fn1[1024],fn2[1024],fn3[1024]
  if(rank<10)strcpy(rankname,"000");
  if(rank<100&&rank>9)strcpy(rankname,"00");
  if(rank<1000&&rank>99)strcpy(rankname,"0");
  if(rank>=1000)rankname[0]='\0';
  sprintf(filename, "%d",rank);
  strcat(rankname,filename);
  nmatin=(uvadd)ndblock/mpisize;
  if(nmatin*mpisize<ndblock)nmatinplus=nmatin+1;
  else nmatinplus=nmatin;
  for(i=0; i<mpisize; i++)if(rank+1<=ndblock-mpisize*nmatin)nmatin++;
  //forint *vnz= (forint *) calloc (nmatin,sizeof(forint));
  //printf("rank %d nmatin %d\n",rank,nmatin);
  begblock[rank]=nmatin;
  for(i=0; i<mpisize; i++) {
    j=rank;
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    j=j2;
    MPI_Bcast(&begblock[j2],1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  j2=begblock[0];
  begblock[0]=0;
  for(i=1; i<mpisize; i++) {
    j1=begblock[i];
    begblock[i]=begblock[i-1]+j2;
    j2=j1;
  }
  //j2=begblock[0];
  //for(i=0; i<mpisize; i++)begblock[i]-=j2;
  //if(rank==1)for(i=0; i<mpisize; i++)printf("begblock %d\n",begblock[i]);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&rowindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&colindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&Cindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&Bindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatA);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatC);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatB);
  CHKERRQ(ierr);
  //PetscInt **ix= (PetscInt**)calloc(nmatin,sizeof(PetscInt*));
  //for (i=0; i<nmatin; i++) ix[i] = (PetscInt*)calloc(ha_ndblocks[i+rank*nmatinplus],sizeof(PetscInt));
  PetscScalar **yi1= (PetscScalar**)calloc(nmatin,sizeof(PetscScalar*));
  for (i=0; i<nmatin; i++) yi1[i] = (PetscScalar*)calloc(ha_ndblocks[i+begblock[rank]],sizeof(PetscScalar));
  //PetscScalar *yi0;
  //Vec *ui1= (Vec*)calloc(nmatin,sizeof(Vec));
  //Mat *submatv= (Mat*)calloc(nmatin,sizeof(Mat));
  //Vec *biui= (Vec*)calloc(nmatin,sizeof(Vec));
  //Mat *submatbivi= (Mat*)calloc(nmatin,sizeof(Mat));
  //for (i=0; i<nmatin; i++) ui[i] = (Vec*)calloc(1,sizeof(Vec));
  //for(i=0; i<nmatin; i++)printf("rank %d proc %d\n",ha_ndblocks[i+rank*nmatinplus],rank);
  sumrowcolin=0;
  for(j=0; j<ndblock; j++)sumrowcolin+=ha_ndblocks[j];
  PetscInt *indicesC= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
  PetscInt *indicesB= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
  PetscScalar *yd = (PetscScalar*)calloc(VecSize-sumrowcolin,sizeof(PetscScalar));
  //for(i=0; i<VecSize-sumrowcolin; i++)yd[i]=-1;
  //PetscScalar *yd0 = (PetscScalar*)calloc(VecSize-sumrowcolin,sizeof(PetscScalar));
  PetscInt *offblock= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
  PetscInt *offblockrow= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
  ha_cgetype *vecbiui= (ha_cgetype *) calloc (VecSize-sumrowcolin,sizeof(ha_cgetype));
  offblock[0]=0;
  for(j=0; j<ndblock; j++)offblock[j+1]=offblock[j]+countvarintra1[j+1]-countvarintra1[j]-ha_ndblocks[j];
  offblockrow[0]=0;
  for(j=0; j<ndblock; j++)offblockrow[j+1]=offblockrow[j]+counteqnoadd[j]-ha_ndblocks[j];
  //for(j=0; j<ndblock; j++)if(rank==0)printf("rank %d row %d eq %d\n",rank,offblockrow[j],counteqnoadd[j]);
  //printf("OK123\n");

  char **fn01= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn01[i] = (char*)calloc(1024,sizeof(char));
  char **fn02= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn02[i] = (char*)calloc(1024,sizeof(char));
  char **fn03= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn03[i] = (char*)calloc(1024,sizeof(char));

  timestr=clock();
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      if(j1<10)strcpy(j1name,"000");
      if(j1<100&&j1>9)strcpy(j1name,"00");
      if(j1<1000&&j1>99)strcpy(j1name,"0");
      if(j1>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j1);
      strcat(j1name,filename);
      strcpy(filename,"_irnv");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn01[j1],filename);
      strcpy(filename,"_keep");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn02[j1],filename);
      strcpy(filename,"_vav");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn03[j1],filename);
    }
  }
  printf("OK012345 rank %d mpis %d\n",rank,mpisize);
  for(j=0; j<mpisize; j++) {
    printf("j %d\n",j);
    for(j3=0; j3<nmatinplus; j3++) {
      if(j3==nmatin)i=j3-1;
      else i=j3;
      //printf("rank %d OK1!!! i %d j %d nmatin %d\n",rank,i,j,nmatin);
      bfirst=counteq[i+begblock[j]];
      bend=ha_ndblocks[i+begblock[j]];
      //if(rank==2)printf("block %d bfirst %d\n",i+begblock[j],bfirst,bend);
      MPI_Bcast(&bfirst,1, MPI_INT,j, PETSC_COMM_WORLD);
      MPI_Bcast(&bend,1, MPI_INT,j, PETSC_COMM_WORLD);
      //if(i<nmatin)bfirst=counteq[i+begblock[j]];
      //else bfirst=VecSize;//counteq[i-1+begblock[j]];
      //if(i<nmatin)bend=ha_ndblocks[i+begblock[j]];
      //else bend=ha_ndblocks[i-1+begblock[j]];
      //bfirst=counteq[i+begblock[j]];
      //bend=ha_ndblocks[i+begblock[j]];
      //if(i<nmatin){
      PetscInt *indices= (PetscInt *) calloc (bend,sizeof(PetscInt));
      ha_cgetype *yi0 = (ha_cgetype*)calloc(bend,sizeof(ha_cgetype));
      //ha_cgetype *yi1 = (ha_cgetype*)calloc(bend,sizeof(ha_cgetype));
      //printf("ha_row 0 %d 588 %d\n",ha_rows[0],ha_rows[588]);
      for(j1=0; j1<bend; j1++) { //ha_ndblocks[i+begblock[j]]
        j2=ha_rows[bfirst+j1]+bfirst;
        if(j2>=Istart&&j2<Iend)indices[j1]=j2;
        else indices[j1]=-1;
        //if(rank==2&&i==21)printf("rank %d Istart %d Iend %d indices %d\n",rank,Istart,Iend,indices[j1]);
      }
      //printf("rank %d OK1!!! i %d j %d\n",rank,i,j);
      VecGetValues(b,bend,indices,yi0);
      //VecView(b,0);
      //for(i=0;i<bend;i++)printf("rank %d yi0 %lf\n",rank,yi0[i]);
      //printf("rank %d OK2!!! i %d j %d nmatin %d\n",rank,i,j,nmatin);

      //if(rank==0)printf("OK rank %d i %d j %d first %d block %d nmatin %d\n",rank,i,j,bfirst,bend,nmatin);
      //MPI_Barrier(PETSC_COMM_WORLD);
      //MPI_Reduce(yi0,yi1,bend, MPI_DOUBLE, MPI_SUM,j,PETSC_COMM_WORLD);
      //forintwrp=bend;
      HaReduce(yi0,(forint)bend,mpisize,rank,j);
      //if(rank==0)printf("OK i %d j %d block %d\n",i,j,bend);
      //MPI_Bcast(yi1,ha_ndblocks[i+begblock[j]]*sizeof(ha_cgetype), MPI_BYTE,0, PETSC_COMM_WORLD);
      //printf("OK1!!! i %d j %d\n",i,j);
      if(rank==j&&i<nmatin) {
        //for(j1=0; j1<bend; j1++) {
        //yi1[i][j1]=yi0[j1];
        //if(rank==0)printf("rank %d i %d indices %d block %d yi %lf\n",rank,i,indices[j1],bend,yi1[i][j1]);
        //}
        memcpy (yi1[i],yi0,bend*sizeof(ha_cgetype));
      }
      free(indices);
      free(yi0);
      //free(yi1);
      //}
      //}
    }
  }
  printf("OK12345 rank %d\n",rank);
  j1=0;
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (ha_ndblocks[i+begblock[rank]],sizeof(PetscInt));
    //printf("ha_row 0 %d 588 %d\n",ha_rows[0],ha_rows[588]);
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_rows[bfirst+j]+bfirst;
      //if(rank==0&&i==1)printf("indices %d rank %d\n",indices[j],ha_ndblocks[i+rank*nmatinplus]);
    }
    //VecGetValues(b,ha_ndblocks[i+rank*nmatinplus],indices,yi[i]);
    //for(j=0;j<ha_ndblocks[i+rank*nmatinplus];j++)printf("rank %d i %d indices %d block %d yi %lf\n",rank,i,indices[j],ha_ndblocks[i+rank*nmatinplus],yi[i][j]);

    bend=counteq[i+begblock[rank]]+counteqnoadd[i+begblock[rank]];
    for(j=offblockrow[i+begblock[rank]]; j<offblockrow[i+1+begblock[rank]]; j++) {
      j1=j-offblockrow[i+begblock[rank]];
      indicesB[j]=ha_rows[bend-j1-1]+bfirst;
      //printf("rank %d rindices %d bend %d j1 %d j %d\n",rank,indicesB[j],ha_rows[bend+j1],j1,j);
    }
    //printf("i %d\n",i);
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
    //printf("i %d\n",i);
    //ISView(rowindices[i],PETSC_VIEWER_STDOUT_WORLD);
    bfirst=countvarintra1[i+begblock[rank]];
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_cols[bfirst+j]+bfirst;
    }
    bend=countvarintra1[i+1+begblock[rank]];
    for(j=offblock[i+begblock[rank]]; j<offblock[i+1+begblock[rank]]; j++) {
      j1=j-offblock[i+begblock[rank]];
      indicesC[j]=ha_cols[bend-j1-1]+bfirst;
      //printf("rank %d cindices %d bend %d j1 %d j %d\n",rank,indicesC[j],ha_cols[bend-j1-1],j1,j);
    }

    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,colindices+i);
    free(indices);
    //ISView(colindices[i],PETSC_VIEWER_STDOUT_WORLD);
    //if(ha_ndblocks[i+rank*nmatinplus]!=0)yi0=realloc (yi0,1*sizeof(PetscScalar));
  }
  printf("0okkkk\n");
  j1=offblock[ndblock];
  for(i=countvarintra1[ndblock]; i<VecSize; i++) {
    indicesC[j1]=i;
    //printf("rank %d indices %d j1 %d\n",rank,indicesC[j1],j1);
    j1++;
  }
  printf("1okkkk\n");
  j1=offblockrow[ndblock];
  for(i=counteq[ndblock-1]+counteqnoadd[ndblock-1]; i<VecSize; i++) {
    indicesB[j1]=i;
    //if(rank==0)printf("rank %d indices %d j1 %d\n",rank,indicesB[j1],j1);
    j1++;
  }
  printf("2okkkk\n");
  //for(i=0; i<VecSize-sumrowcolin; i++)if(rank==0)printf("0rank %d bindx %d\n",rank,indicesB[i]);
  for(i=0; i<mpisize; i++) {
    j=offblockrow[begblock[rank]];
    j1=offblockrow[nmatin+begblock[rank]]-offblockrow[begblock[rank]];
    //printf("rank %d i %d j %d j1 %d\n",rank,i,j,j1);
    MPI_Bcast(&j,1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&j1,1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&indicesB[j],j1, MPI_INT,i, PETSC_COMM_WORLD);
    //if(rank==i)VecGetValues(b,j1,&indicesB[j],&yd[j]);
  }
  //printf("rank0 %d\n",rank);
  //for(i=0; i<offblockrow[ndblock]; i++)if(rank==1)printf("1rank %d bindx %d\n",rank,indicesB[i]);
  //for(i=0; i<VecSize-sumrowcolin; i++)if(rank==2)printf("indiceB %d\n",indicesB[i]);
  printf("rank1 %d\n",rank);
  //for(i=0; i<VecSize-sumrowcolin; i++){if(rank==1)printf("i %d\n",i);VecGetValues(b,1,&indicesB[i],&yd[i]);}
  //for(i=0; i<mpisize; i++)if(rank==i)

  //if(rank==mpisize-1){
  //bend=VecSize-counteq[ndblock-1]-counteqnoadd[ndblock-1];
  //VecGetValues(b,bend,&indicesB[j1+j],&yd[j1+j]);
  //}

  for(i=0; i<VecSize-sumrowcolin; i++)if(indicesB[i]>=Istart) {
      j=i;
      break;
    }
  for(j2=i; j2<VecSize-sumrowcolin; j2++)if(indicesB[j2]>=Iend) {
      j1=j2;
      break;
    }
  if(j2==VecSize-sumrowcolin)j1=VecSize-sumrowcolin;
  printf("rank %d j %d j1 %d istart %d iend %d\n",rank,j,j1,Istart,Iend);
  VecGetValues(b,j1-j,&indicesB[j],&yd[j]);
  //MPI_Barrier(PETSC_COMM_WORLD);
  //if(rank==mpisize-1)VecGetValues(b,VecSize-sumrowcolin-bend,indicesB+bend,yd+bend);
  printf("rank2 %d\n",rank);
  //MPI_Reduce(yd,yd0,VecSize-sumrowcolin, MPI_DOUBLE, MPI_SUM,0,PETSC_COMM_WORLD);
  forintwrp=VecSize-sumrowcolin;
  HaReduce(yd,forintwrp,mpisize,rank,mpisize-1);
  if(rank!=mpisize-1)free(yd);
  //if(rank!=0)free(yd0);
  //printf("Vec %d sumrow %d\n",VecSize,sumrowcolin);
  //for(i=0; i<VecSize-sumrowcolin; i++)if(rank==0)printf("rank %d yd %f\n",rank,yd0[i]);
  for(i=0; i<mpisize; i++) {
    j=offblock[begblock[rank]];
    j1=offblock[nmatin+begblock[rank]]-offblock[begblock[rank]];
    //printf("rank %d i %d j %d j1 %d\n",rank,i,j,j1);
    MPI_Bcast(&j,1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&j1,1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&indicesC[j],j1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  printf("OK132\n");
  //for(i=0; i<mpisize; i++) {
  //j=offblockrow[rank*nmatinplus];
  //j1=offblockrow[nmatin+rank*nmatinplus]-offblockrow[rank*nmatinplus];
  //printf("rank %d i %d j %d j1 %d\n",rank,i,j,j1);
  //MPI_Bcast(&j,sizeof(PetscInt), MPI_BYTE,i, PETSC_COMM_WORLD);
  //MPI_Bcast(&j1,sizeof(PetscInt), MPI_BYTE,i, PETSC_COMM_WORLD);
  //MPI_Bcast(&indicesB[j],j1*sizeof(PetscInt), MPI_BYTE,i, PETSC_COMM_WORLD);
  //}
  for(i=0; i<nmatin; i++) {
    ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesC,PETSC_COPY_VALUES,Cindices+i);
    ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesB,PETSC_COPY_VALUES,Bindices+i);
  }
  free(indicesB);
  //for(i=0; i<VecSize-sumrowcolin; i++)printf("rank %d indiceB %d\n",rank,indicesB[i]);
  //for(i=0; i<VecSize-sumrowcolin; i++)printf("rank %d indiceC %d\n",rank,indicesC[i]);
  free(offblock);
  free(offblockrow);
  printf("Begin partitioning rank %d\n",rank);
  //if(rank==0)ISView(Bindices[0],PETSC_VIEWER_STDOUT_SELF);
  //if(rank==0)ISView(Cindices[0],PETSC_VIEWER_STDOUT_SELF);
  //if(rank==0)ISView(rowindices[0],PETSC_VIEWER_STDOUT_SELF);
  //if(rank==0)ISView(colindices[0],PETSC_VIEWER_STDOUT_SELF);
  ierr = MatCreateSubMatrices(A,nmatin,rowindices,colindices,MAT_INITIAL_MATRIX,&submatA);
  printf("A\n");
  ierr = MatCreateSubMatrices(A,nmatin,rowindices,Cindices,MAT_INITIAL_MATRIX,&submatC);
  printf("C\n");
  ierr = MatCreateSubMatrices(A,nmatin,Bindices,colindices,MAT_INITIAL_MATRIX,&submatB);
  printf("B\n");
  /*int *nzbs;//(uvadd *) calloc (1,sizeof(int));
  int *bns;//=(int *) calloc (1,sizeof(int));
  int *bms;//=(int *) calloc (1,sizeof(int));
  int *irnbs;//=(int *) calloc (1,sizeof(int));
  int *jcnbs;//=(int *) calloc (1,sizeof(int));
  PetscScalar *valbs;//=(PetscScalar *) calloc (1,sizeof(PetscScalar));

  IS *colindices1,*rowindices1;
  int nz1,indxbsize=VecSize-sumrowcolin;
  for(i=0; i<mpisize-1; i++) {
    j0=nmatin;
    MPI_Bcast(&j0,sizeof(j0), MPI_BYTE,i, PETSC_COMM_WORLD);
    //printf("i0 %d rank %d\n",i,rank);
    if(rank==mpisize-1) {
      //IS *colindices= (IS *) calloc (j0,sizeof(**IS));
      PetscMalloc(j0*sizeof(IS **),&colindices1);
      PetscMalloc(j0*sizeof(IS **),&rowindices1);
      for(j=0; j<j0; j++) {
        PetscInt *indices= (PetscInt *) calloc (ha_ndblocks[j+begblock[i]],sizeof(PetscInt));
        bfirst=countvarintra1[j+begblock[i]];
        for(j1=0; j1<ha_ndblocks[j+begblock[i]]; j1++) {
          indices[j1]=ha_cols[bfirst+j1]+bfirst;
          //printf("indice %d j1 %d size %d\n",indices[j1],j1,ha_ndblocks[j+begblock[i]]);
        }
        ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[j+begblock[i]],indices,PETSC_COPY_VALUES,colindices1+j);
        //ISView(rowindices[j],PETSC_VIEWER_STDOUT_SELF);
        free(indices);
        ISCreateGeneral(PETSC_COMM_SELF,indxbsize,indicesB,PETSC_COPY_VALUES,rowindices1+j);
      }
      PetscMalloc(j0*sizeof(Mat *),&submatB1);
      //printf("j0 %d\n",j0);
    }
    j1=j0;
    if(rank!=mpisize-1)j1=0;
    MatCreateSubMatrices(A,j1,rowindices1,colindices1,MAT_INITIAL_MATRIX,&submatB1);
    if(rank==mpisize-1) {
      for(j=0; j<j0; j++) {
        ISDestroy(colindices1+j);
        ISDestroy(rowindices1+j);
      }
      PetscFree(colindices1);
      PetscFree(rowindices1);
      //MatView(submatB[0],0);
      nz=0;
      nz1=0;
      for(j=0; j<j0; j++) {
        Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatB1[j]->data;
        //nzbs[j]=ab->nz;
        nz+=ab->nz;//nzbs[j];
        nz1+=submatB1[j]->cmap->n;//nzbs[j];
      }
      if(i<mpisize-1) {
        MPI_Send(&nz,1, MPI_INT, i, 0, PETSC_COMM_WORLD);
        MPI_Send(&nz1,1, MPI_INT, i, 7, PETSC_COMM_WORLD);
      }
      //printf("i %d rank %d nz %d\n",i,rank,nz);
    }
    if(rank==i) {
      MPI_Recv(&nz,1, MPI_INT, mpisize-1, 0, PETSC_COMM_WORLD,&status);
      MPI_Recv(&nz1,1, MPI_INT, mpisize-1, 7, PETSC_COMM_WORLD,&status);
    }
    if(rank==i||rank==mpisize-1) {
      //*nzbs= calloc (*nzbs,j0*sizeof(int));
      //*bns= calloc (*bns,j0*sizeof(int));
      //*bms= calloc (*bms,j0*sizeof(int));
      //*jcnbs= calloc (*jcnbs,nz*sizeof(int));
      //*irnbs= calloc (*irnbs,nz1*sizeof(int));
      //*valbs= calloc (*valbs,nz*sizeof(PetscScalar));
      nzbs= (int*)calloc (j0,sizeof(int));
      bns= (int*)calloc (j0,sizeof(int));
      bms= (int*)calloc (j0,sizeof(int));
      jcnbs= (int*)calloc (nz,sizeof(int));
      irnbs= (int*)calloc (nz1,sizeof(int));
      valbs= (PetscScalar*)calloc (nz,sizeof(PetscScalar));
      printf("i5 %d rank %d nz %d j0 %d nz1 %d\n",i,rank,nz,j0,nz1);
    }
    if(rank==mpisize-1) {
      j3=0;
      j4=0;
      for(j=0; j<j0; j++) {
        MatTranspose(submatB1[j],MAT_INITIAL_MATRIX,&submatBT);
        MatDestroy_SeqAIJ(submatB1[j]);
        Mat_SeqAIJ         *abt=(Mat_SeqAIJ*)submatBT->data;
        ai= abt->i;
        aj= abt->j;
        vals=abt->a;
        nzbs[j]=abt->nz;
        bns[j]=submatBT->rmap->n;
        bms[j]=submatBT->cmap->n;
        //printf("j %d btn %d bms %d nz %d\n",j,submatBT->cmap->n,(*bms)[j],(*nzbs)[j]);
        for(j1=0; j1<bns[j]; j1++) {
          irnbs[j4]=ai[j1];
          j4++;
        }
        for(j1=0; j1<nzbs[j]; j1++) {
          jcnbs[j3]=aj[j1];
          //printf("jcn %d\n",(*jcnbs)[j3]);
          valbs[j3]=vals[j1];
          j3++;
        }
        MatDestroy_SeqAIJ(submatBT);
      }
      free(submatB1);
    }
    if(i<mpisize-1) {
      if(rank==mpisize-1) {
        MPI_Send(valbs,nz, MPI_DOUBLE,i, 31, PETSC_COMM_WORLD);
        MPI_Send(irnbs,nz1, MPI_INT,i, 32, PETSC_COMM_WORLD);
        MPI_Send(jcnbs,nz, MPI_INT,i, 33, PETSC_COMM_WORLD);
        MPI_Send(bns,j0, MPI_INT,i, 34, PETSC_COMM_WORLD);
        MPI_Send(bms,j0, MPI_INT,i, 35, PETSC_COMM_WORLD);
        MPI_Send(nzbs,j0, MPI_INT,i, 36, PETSC_COMM_WORLD);
        for(j1=0; j1<j0; j1++)printf("nzzzzzzsend %d rank %d i %d j0 %d j1 %d\n",nzbs[j1],rank,i,j0,j1);
      }
      if(rank==i) {
        MPI_Recv(valbs,nz, MPI_DOUBLE,mpisize-1, 31, PETSC_COMM_WORLD,&status);
        MPI_Recv(irnbs,nz1, MPI_INT,mpisize-1, 32, PETSC_COMM_WORLD,&status);
        MPI_Recv(jcnbs,nz, MPI_INT,mpisize-1, 33, PETSC_COMM_WORLD,&status);
        MPI_Recv(bns,j0, MPI_INT,mpisize-1, 34, PETSC_COMM_WORLD,&status);
        MPI_Recv(bms,j0, MPI_INT,mpisize-1, 35, PETSC_COMM_WORLD,&status);
        MPI_Recv(nzbs,j0, MPI_INT,mpisize-1, 36, PETSC_COMM_WORLD,&status);
        //for(j1=0; j1<j0; j1++)printf("nzzzzzz %d rank %d\n",(*nzbs)[j1],rank);
        if(rank==0)for(j1=0; j1<j0; j1++)printf("nzzzzzzrecv %d rank %d i %d j0 %d j1 %d\n",nzbs[j1],rank,i,j0,j1);
      }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    //if(rank==0)for(j1=0; j1<j0; j1++)printf("nzzzzzz %d rank %d i %d\n",nzbs[j1],rank,i);
  }

  i=mpisize-1;
  j0=nmatin;
  if(rank==mpisize-1) {
    //IS *colindices= (IS *) calloc (j0,sizeof(**IS));
    PetscMalloc(j0*sizeof(IS **),&colindices1);
    PetscMalloc(j0*sizeof(IS **),&rowindices1);
    for(j=0; j<j0; j++) {
      PetscInt *indices= (PetscInt *) calloc (ha_ndblocks[j+begblock[i]],sizeof(PetscInt));
      bfirst=countvarintra1[j+begblock[i]];
      for(j1=0; j1<ha_ndblocks[j+begblock[i]]; j1++) {
        indices[j1]=ha_cols[bfirst+j1]+bfirst;
        //printf("indice %d j1 %d size %d\n",indices[j1],j1,ha_ndblocks[j+begblock[i]]);
      }
      ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[j+begblock[i]],indices,PETSC_COPY_VALUES,colindices1+j);
      //ISView(rowindices[j],PETSC_VIEWER_STDOUT_SELF);
      free(indices);
      ISCreateGeneral(PETSC_COMM_SELF,indxbsize,indicesB,PETSC_COPY_VALUES,rowindices1+j);
    }
    PetscMalloc(j0*sizeof(Mat *),&submatB1);
    //printf("j0 %d\n",j0);
  }
  j1=j0;
  if(rank!=mpisize-1)j1=0;
  MatCreateSubMatrices(A,j1,rowindices1,colindices1,MAT_INITIAL_MATRIX,&submatB1);
  if(rank==mpisize-1) {
    for(j=0; j<j0; j++) {
      ISDestroy(colindices1+j);
      ISDestroy(rowindices1+j);
    }
    PetscFree(colindices1);
    PetscFree(rowindices1);
    nz=0;
    nz1=0;
    for(j=0; j<j0; j++) {
      Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatB1[j]->data;
      //nzbs[j]=ab->nz;
      nz+=ab->nz;//nzbs[j];
      nz1+=submatB1[j]->cmap->n;//nzbs[j];
    }
    //*nzbs= realloc (*nzbs,j0*sizeof(int));
    //*bns= realloc (*bns,j0*sizeof(int));
    //*bms= realloc (*bms,j0*sizeof(int));
    //*jcnbs= realloc (*jcnbs,nz*sizeof(int));
    //*irnbs= realloc (*irnbs,nz1*sizeof(int));
    //*valbs= realloc (*valbs,nz*sizeof(PetscScalar));
    nzbs= (int*)calloc (j0,sizeof(int));
    bns= (int*)calloc (j0,sizeof(int));
    bms= (int*)calloc (j0,sizeof(int));
    jcnbs= (int*)calloc (nz,sizeof(int));
    irnbs= (int*)calloc (nz1,sizeof(int));
    valbs= (PetscScalar*)calloc (nz,sizeof(PetscScalar));
    j3=0;
    j4=0;
    for(j=0; j<j0; j++) {
      MatTranspose(submatB1[j],MAT_INITIAL_MATRIX,&submatBT);
      MatDestroy_SeqAIJ(submatB1[j]);
      Mat_SeqAIJ         *abt=(Mat_SeqAIJ*)submatBT->data;
      ai= abt->i;
      aj= abt->j;
      vals=abt->a;
      nzbs[j]=abt->nz;
      bns[j]=submatBT->rmap->n;
      bms[j]=submatBT->cmap->n;
      //printf("j %d btn %d bms %d nz %d\n",j,submatBT->cmap->n,(*bms)[j],(*nzbs)[j]);
      for(j1=0; j1<bns[j]; j1++) {
        irnbs[j4]=ai[j1];
        j4++;
      }
      for(j1=0; j1<nzbs[j]; j1++) {
        jcnbs[j3]=aj[j1];
        //printf("jcn %d\n",(*jcnbs)[j3]);
        valbs[j3]=vals[j1];
        j3++;
      }
      MatDestroy_SeqAIJ(submatBT);
    }
    free(submatB1);
  }*/


  //HaMatGetSubMats(A,VecSize-sumrowcolin,indicesB,ha_ndblocks,begblock,countvarintra1,ha_cols,mpisize,rank,nmatin,mpisize-1,&nzbs,&bns,&bms,&irnbs,&jcnbs,&valbs);
  //if(rank==0)for(i=0;i<1159;i++)printf("i %d jcn %d val %lf\n",i,jcnbs[i],valbs[i]);
  //if(rank==2)for(i=0; i<nmatin; i++)printf("i12345 %d nz %d\n",i,nzbs[i]);
  //printf("rank %d nzbs %d\n",rank,bms[0]);
  //Mat *submatB = (Mat*)calloc(nmatin,sizeof(Mat));
  //for(i=0; i<nmatin; i++)MatTranspose(submatB1[i],MAT_INITIAL_MATRIX,&submatB[i]);
  if(rank==mpisize-1)proc1=1;
  ierr = PetscMalloc(proc1*sizeof(Mat *),&submatD);
  CHKERRQ(ierr);
  printf("D\n");
  MatCreateSubMatrices(A,proc1,Bindices,Cindices,MAT_INITIAL_MATRIX,&submatD);
  printf("End partitionning rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  //CHKERRQ(ierr);
  //PetscViewerASCIIOpen(PETSC_COMM_SELF, "A.txt", &viewer);
  //PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_DENSE );
  //ierr = MatView(A,viewer);
  //CHKERRQ(ierr);
  //PetscViewerDestroy(&viewer);
  ierr = MatDestroy(&A);
  CHKERRQ(ierr);
  ierr = VecDestroy(&b);
  CHKERRQ(ierr);
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&rowindices[i]);
    CHKERRQ(ierr);
    //ierr = ISDestroy(colindices+i);
    //CHKERRQ(ierr);
    ierr = ISDestroy(&Cindices[i]);
    CHKERRQ(ierr);
    ierr = ISDestroy(&Bindices[i]);
    CHKERRQ(ierr);
  }
  PetscFree(rowindices);
  PetscFree(Cindices);
  PetscFree(Bindices);
  //PetscViewerASCIIOpen(PETSC_COMM_SELF, "b.txt", &viewer);
  //ierr = VecView(b,0);
  //CHKERRQ(ierr);
  //PetscViewerDestroy(&viewer);
  if(rank<0) {
    //PetscViewerASCIIOpen(PETSC_COMM_SELF, "Asub.txt", &viewer);
    //PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_DENSE );
    //ierr = MatView(submatA[0],viewer);
    //CHKERRQ(ierr);
    //PetscViewerDestroy(&viewer);
    PetscViewerASCIIOpen(PETSC_COMM_SELF, "C.txt", &viewer);
    //PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_DENSE );
    ierr = MatView(submatC[0],viewer);
    //CHKERRQ(ierr);
    PetscViewerDestroy(&viewer);
    PetscViewerASCIIOpen(PETSC_COMM_SELF, "B.txt", &viewer);
    //PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_DENSE );
    ierr = MatView(submatB[0],viewer);
    //CHKERRQ(ierr);
    PetscViewerDestroy(&viewer);
    PetscViewerASCIIOpen(PETSC_COMM_SELF, "D.txt", &viewer);
    //PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_DENSE );
    ierr = MatView(submatD[0],viewer);
    //CHKERRQ(ierr);
    PetscViewerDestroy(&viewer);
  }
  printf("OK1 rank %d\n",rank);
  //free(ix);
  //printf("rank0 %d vec %lf\n",rank,vecbivi[22800]);
  //ha_cgetype *vecbivi= (ha_cgetype *) calloc ((VecSize-sumrowcolin)*(VecSize-sumrowcolin),sizeof(ha_cgetype));
  ha_cgetype *xi1 = (ha_cgetype*)calloc(sumrowcolin,sizeof(ha_cgetype));
  int insizes=17;
  int *insize=(int *) calloc (insizes*nmatin,sizeof(int));
  int *bivinzrow=(int *) calloc (VecSize-sumrowcolin,sizeof(int));
  int *bivinzcol=(int *) calloc (VecSize-sumrowcolin,sizeof(int));
  //int bivinz0=0;
  //if(rank=mpisize-1)bivinz0=VecSize-sumrowcolin;
  //int *bivinzrow0=(int *) calloc (VecSize-sumrowcolin,sizeof(int));
  //int *bivinzcol0=(int *) calloc (VecSize-sumrowcolin,sizeof(int));
  //j1nz=0;
  //j1irnbs=0;
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      MatTranspose(submatC[j1],MAT_INITIAL_MATRIX,&submatBT);
      Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatBT->data;//*aa=subA->data;
      ai= ac->i;
      nz=ac->nz;
      ncol=submatC[j1]->cmap->n;
      for(i=0; i<ncol-1; i++)if(ai[i]!=ai[i+1])bivinzcol[i]++;
      if(ai[i]<nz)bivinzcol[i]++;
      MatDestroy_SeqAIJ(submatBT);
      Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatB[j1]->data;//*aa=subA->data;
      ai= ab->i;
      nz=ab->nz;
      nrow=submatB[j1]->rmap->n;
      for(i=0; i<nrow-1; i++)if(ai[i]!=ai[i+1])bivinzrow[i]++;
      if(ai[i]<nz)bivinzrow[i]++;
    }
  }
  int bivirowsize=1,bivicolsize=1;
  forint vecbivisize=0;
  //MPI_Allreduce(bivinzrow,bivinzrow0,VecSize-sumrowcolin, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
  //MPI_Allreduce(bivinzcol,bivinzcol0,VecSize-sumrowcolin, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
  //if(rank==mpisize-1){
  //memcpy(bivinzcol,bivinzcol0,(VecSize-sumrowcolin)*sizeof(int));
  //memcpy(bivinzrow,bivinzrow0,(VecSize-sumrowcolin)*sizeof(int));
  //}
  for(i=0; i<VecSize-sumrowcolin; i++) {
    if(bivinzrow[i]>0) {
      bivinzrow[i]=bivirowsize;
      bivirowsize++;
    }
    if(bivinzcol[i]>0) {
      bivinzcol[i]=bivicolsize;
      bivicolsize++;
    }
  }
  bivirowsize--;
  bivicolsize--;
  int *bivinzrow1=(int *) calloc (bivirowsize,sizeof(int));
  int *bivinzcol1=(int *) calloc (bivicolsize,sizeof(int));
  j=0;
  j2=0;
  for(i=0; i<VecSize-sumrowcolin; i++) {
    if(bivinzrow[i]>0) {
      bivinzrow1[j]=i;
      j++;
    }
    if(bivinzcol[i]>0) {
      bivinzcol1[j2]=i;
      j2++;
    }
  }
  //if(rank==mpisize-1){
  //memcpy(bivinzcol0,bivinzcol,(VecSize-sumrowcolin)*sizeof(int));
  //memcpy(bivinzrow0,bivinzrow,(VecSize-sumrowcolin)*sizeof(int));
  //}
  //MPI_Bcast(bivinzrow0,(VecSize-sumrowcolin)*sizeof(int), MPI_BYTE,mpisize-1, PETSC_COMM_WORLD);
  //MPI_Bcast(bivinzcol0,(VecSize-sumrowcolin)*sizeof(int), MPI_BYTE,mpisize-1, PETSC_COMM_WORLD);
  //for(i=0;i<VecSize-sumrowcolin;i++)printf("rank %d bivirow %d 0 %d\n",rank,bivinzrow[i],bivinzrow0[i]);
  //printf("f %d c %d r %d\n",VecSize-sumrowcolin,bivirowsize,bivicolsize);
  //if(rank==0)for(i=0;i<VecSize-sumrowcolin;i++)printf("bivirow %d\n",bivinzcol0[i]);
  vecbivisize=bivirowsize*bivicolsize;
  ha_cgetype *vecbivi= (ha_cgetype *) calloc (vecbivisize,sizeof(ha_cgetype));
  printf("vecbivisize %ld rank %d\n",vecbivisize,rank);

  ha_cgetype *xi1point;
  uvadd xi1indx=0;
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      time(&timestr);//=clock();
      printf("OK1 j1 %d rank %d\n",j1,rank);
      Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatA[j1]->data;//*aa=subA->data;
      //printf("OK1 j1 %d\n",j1);
      MatTranspose(submatC[j1],MAT_INITIAL_MATRIX,&submatBT);
      //MatDestroy_SeqAIJ(submatC[j1]);
      Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatBT->data;//*aa=subA->data;
      //printf("0 rank %d j1 %d vec %lf\n",rank,j1,vecbivi[22800]);
      ai= aa->i;
      aj= aa->j;
      vals=aa->a;
      nz=aa->nz;
      nrow=submatA[j1]->rmap->n;
      ncol=submatA[j1]->cmap->n;
      //printf("nz %d\n",nz);
      //la=2;
      lasize=ceil((laA/100.0)*nz);
      int *irn=(int *) calloc (lasize,sizeof(int));
      int *jcn=(int *) calloc (lasize,sizeof(int));
      ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
      //ha_cgetype *xi = (ha_cgetype*)calloc(ha_ndblocks[j1+begblock[rank]],sizeof(ha_cgetype));
      //int *insize=(int *) calloc (12,sizeof(int));
      //long long int *insize=(long long int *) calloc (10,sizeof(long long int));
      insize[j1*insizes+13]=bivirowsize;
      insize[j1*insizes+14]=bivicolsize;
      insize[j1*insizes]=nrow;
      insize[j1*insizes+1]=ncol;
      insize[j1*insizes+2]=nz;
      aic= ac->i;
      ajc= ac->j;
      valsc=ac->a;
      nzc=ac->nz;
      nrowc=submatBT->rmap->n;
      ncolc=submatBT->cmap->n;
      printf("nrowc %d ncolc %d nzc %d\n",nrowc,ncolc,nzc);
      int *irnc=(int *) calloc (nrowc,sizeof(int));
      //int *irnc1=(int *) calloc (nrowc,sizeof(int));
      int *jcnc=(int *) calloc (nzc,sizeof(int));
      ha_cgetype *valuesc= (ha_cgetype *) calloc (nzc,sizeof(ha_cgetype));
      insize[j1*insizes+3]=nrowc;
      insize[j1*insizes+4]=ncolc;
      //for(i=0; i<nrow; i++)if(rank==1)printf("0aici %d\n",ai[i]);
      if(nrow==1)for(i=0; i<nz; i++)irn[i]=1;
      else {
        if(ai[1]>0)for(j=ai[0]; j<ai[1]; j++)irn[j]=1;
        for(i=1; i<nrow-1; i++)if(ai[i+1]>ai[i]) {
            for(j=ai[i]; j<ai[i+1]; j++) {
              irn[j]=i+1;
              //if(rank==1&&(i==5||i==6))printf("aici %d aicip1 %d j %d irn %d\n",aic[i],aic[i+1],j,irn[j]);
            }
          }
      }
      if(ai[nrow-1]<nz)for(j=ai[nrow-1]; j<nz; j++) {
          irn[j]=nrow;
        }
      //j=0;
      for(i=0; i<nz; i++) {
        //if(vals[i]!=0) {
        //irn[i]=irn[i];
        jcn[i]=aj[i]+1;
        //values[i]=vals[i];
        //j++;
        //}
        //printf("irn %d jcn %d val %lf\n",irn[i],jcn[i],vals[i]);
      }
      memcpy (values,vals,nz*sizeof(ha_cgetype));
      printf("nz %d j %d\n",nz,j);
      //nz=j;
      //jcn=realloc(jcn,la*nz*sizeof(int));
      //irn=realloc(irn,la*nz*sizeof(int));
      //values=realloc(values,la*nz*sizeof(ha_cgetype));
      insize[j1*insizes+2]=nz;
      /*for(i=0; i<nz; i++) {
        jcn[i]=i+1;
        values[i]=vals[i];
        irn[i]=i+1;
        //if(rank==1)printf("irn %d jcn %d val %lf\n",irn[i],jcn[i],values[i]);
      }*/
      /*for(i=0; i<nrow-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
          irn[j]=i+1;
        }
      for(j=ai[nrow-1]; j<nz; j++) {
        irn[j]=nrow;
      }
      j=0;
      for(i=0; i<nz; i++) {
        jcn[i]=aj[i]+1;
        values[i]=vals[i];
      }*/
      //for(i=0; i<nrowc-1; i++)if(rank==1)printf("0aici %d\n",irnc[i]);
      for(i=0; i<nrowc-1; i++)irnc[i]=aic[i+1]-aic[i];//printf("i %d irnc %d\n",i,irnc[i]);}
      irnc[i]=nzc-aic[i];//printf("i %d irnc %d\n",i,irnc[i]);
      for(i=0; i<nzc; i++) {
        jcnc[i]=ajc[i]+1;
        //printf("jnc %d valsc %lf\n",jcnc[i],valsc[i]);
      }
      memcpy (valuesc,valsc,nzc*sizeof(ha_cgetype));
      //memcpy (jcnc,ajc,nzc*sizeof(int));
      //for(i=0; i<nzc; i++)jcnc++;
      /*if(nrowc==1)for(i=0; i<nzc; i++)irnc[i]=1;
      else {
        if(aic[1]>0)for(j=aic[0]; j<aic[1]; j++)irnc[j]=1;
        for(i=1; i<nrowc-1; i++)if(aic[i+1]>aic[i]) {
            for(j=aic[i]; j<aic[i+1]; j++) {
              //printf("aici %d aicip1 %d j %d\n",aic[i],aic[i+1],j);
              irnc[j]=i+1;
            }
          }
      }
      if(aic[nrowc-1]<nzc)for(j=aic[nrowc-1]; j<nzc; j++) {
          irnc[j]=nrowc;
        }*/
      //j=0;
      //j3=0;
      //for(i=0; i<nrowc; i++) {
      //j4=0;
      //for(j2=0; j2<irnc[i]; j2++) {
      //j5=j2+j3;
      //if(valsc[j5]!=0) {
      //irnc[j]=irnc[j2+j3];
      //jcnc[j]=ajc[j5]+1;
      //printf("jcnc %d\n",jcnc[j]);
      //valuesc[j]=valsc[j5];
      //j++;
      //} else j4+=1;
      //}
      //j3+=irnc[i];
      //irnc[i]-=j4;
      //}
      //nzc=j;
      //jcnc=realloc(jcnc,nzc*sizeof(int));
      //irnc=realloc(irnc,nzc*sizeof(int));
      //valuesc=realloc(valuesc,nzc*sizeof(ha_cgetype));
      insize[j1*insizes+5]=nzc;
      MatDestroy(&submatA[j1]);
      MatDestroy_SeqAIJ(submatBT);
      //for(i=0; i<nzc; i++)printf("irnc %d\n",irnc[i]);
      //for(i=0; i<10; i++) printf("x %f irn %d jcn %d n %d m %d\n",values[i],irn[i],jcn[i],nrow,ncol);
      //extern void spec51_rank_(forint *INSIZE,forint *IRN, forint *JCN, ha_cgetype *VA);
      //nzv=nrowc/7;//(VecSize-sumrowcolin)*ha_ndblocks[j1+rank*nmatinplus];
      //nzv=nzv*ncolc;
      //insize[6]=nrowc;
      //insize[7]=ncolc;
      //insize[8]=nzv;
      insize[j1*insizes+9]=laA;
      insize[j1*insizes+10]=rank;
      insize[j1*insizes+11]=j1;
      //printf("nzc %d OK????\n",nzc);
      //forint *irnv=(forint *) calloc (nzv,sizeof(forint));
      //forint *jcnv=(forint *) calloc (nzv,sizeof(forint));
      //ha_cgetype *valuesv= (ha_cgetype *) calloc (nzv,sizeof(ha_cgetype));
      MatGetSize(submatB[j1],&nrowb,&ncolb);
      //MatGetSize(submatB[j1],&nzv,&nnzmax);
      printf("rank1 %d nzv size %d\n",rank,nzv);
      MatTranspose(submatB[j1],MAT_INITIAL_MATRIX,&submatBT);
      MatDestroy(&submatB[j1]);
      Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatBT->data;//*aa=subA->data;
      ai= ab->i;//irnbs+j1irnbs;
      aj= ab->j;//jcnbs+j1nz;
      vals=ab->a;//valbs+j1nz;
      nz=ab->nz;//nzbs[j1];
      insize[j1*insizes+6]=ncolb;
      insize[j1*insizes+7]=nrowb;
      insize[j1*insizes+8]=nz;
      xi1point=xi1+xi1indx;
      //for(i=0;i<10;i++)printf("vecbui0 %lf\n",vecbivi[i]);//bivirowsize*bivicolsize
      spec48_msol_(insize+j1*insizes,irn,jcn,values,yi1[j1],xi1point,irnc,jcnc,valuesc,ai,aj,vals,vecbivi,bivinzrow,bivinzcol);
      //for(i=0;i<10;i++)printf("vecbui1 %lf\n",vecbivi[i]);//
      printf("rank %d nzv1 %d OK????\n",rank,nzv);
      //if(rank==0)for(i=0;i<ncolb;i++)printf("rank %d i %d xi %lf\n",rank,i,xi1point[i]);
      //vnz[j1]=insize[8];
      free(yi1[j1]);
      free(irnc);
      free(jcnc);
      free(valuesc);
      //if(insize[8]>nzv)printf("Errors: Far too many nonzeros. Please increase nz ratio!");
      //if(rank==1)spec48_ssol_(insize,irn,jcn,values,yi[j1],xi);
      /*MatDuplicate(submatA[j1],MAT_DO_NOT_COPY_VALUES,&submatAA);
      for(i=0; i<insize[2]; i++) {
        irn[i]=irn[i]-1;
        jcn[i]=jcn[i]-1;
        printf("i %d r %d c %d v %lf\n",i,irn[i],jcn[i],values[i]);
        MatSetValues(submatAA,1,(PetscInt)&irn[i],1,(PetscInt)&jcn[i],&values[i],INSERT_VALUES);
      }
      if(rank==0){
      MatAssemblyBegin(submatAA,MAT_FINAL_ASSEMBLY); MatAssemblyEnd(submatAA,MAT_FINAL_ASSEMBLY);
      PetscViewerASCIIOpen(PETSC_COMM_SELF, "AsubAA.txt", &viewer);
      //PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_DENSE );
      ierr = MatView(submatAA,viewer);
      CHKERRQ(ierr);
      PetscViewerDestroy(&viewer);
      }
      //MatDestroy_SeqAIJ(submatAA);*/
      //spec48_ssol_(insize,irn,jcn,values,yi[j1],xi);
      //printf("sizeof forint %d OK!!!!\n",sizeof(forint));
      free(irn);
      free(jcn);
      free(values);
      //Multiply Bi by ui:
      for(i=0; i<ncolb-1; i++) {
        for(j=ai[i]; j<ai[i+1]; j++) {
          vecbiui[aj[j]]+=vals[j]*xi1point[i];//xi[i];
        }
      }
      for(j=ai[i]; j<nz; j++) {
        vecbiui[aj[j]]+=vals[j]*xi1point[i];//xi[i];
      }
      //if(rank==0)for(i=0;i<ncolb;i++)printf("rank %d i %d xi %lf\n",rank,i,vecbiui[i]);
      //for(i=0; i<ha_ndblocks[j1+rank*nmatinplus]; i++)if(xi[i]>=0)printf("i %d xi %lf\n",i,xi[i]);
      //VecCreateSeq(PETSC_COMM_SELF,ha_ndblocks[j1+begblock[rank]],&ui1[j1]);
      //VecGetArray(ui1[j1],&vals);
      //for(i=0; i<ha_ndblocks[j1+begblock[rank]]; i++)if(xi[i]!=0)vals[i]=xi[i];//VecSetValue(ui1[j1],i,xi[i], INSERT_VALUES);
      //VecRestoreArray(ui1[j1],&vals);
      //VecCreateSeqWithArray(PETSC_COMM_SELF,1,ha_ndblocks[j1+begblock[rank]],xi,&ui1[j1]);

      //VecSetValues(ui1[j1],ha_ndblocks[j1+begblock[rank]],xi, INSERT_VALUES);
      //i=0;
      //VecSetBlockSize(ui1[j1],ha_ndblocks[j1+begblock[rank]]);
      //VecSetValuesBlocked(ui1[j1],1,&i,xi,INSERT_VALUES);

      //free(xi);
      //VecView(ui[j1],0);
      //VecAssemblyBegin(ui1[j1]);
      //VecAssemblyEnd(ui1[j1]);
      //PetscInt *nnz=(PetscInt *) calloc (ha_ndblocks[j1+rank*nmatinplus],sizeof(PetscInt));
      //PetscInt *nnzc=(PetscInt *) calloc ((VecSize-sumrowcolin),sizeof(PetscInt));
      //printf("inside %d\n",insize[8]);
      //printf("rank %d nzv size %d\n",rank,nzv);
      //if(j1<10)strcpy(j1name,"000");
      //if(j1<100&&j1>9)strcpy(j1name,"00");
      //if(j1<1000&&j1>99)strcpy(j1name,"0");
      //sprintf(filename, "%d",j1);
      //strcat(j1name,filename);
      //strcpy(filename,"_irnv");
      //strcat(filename,rankname);
      //strcat(filename,j1name);
      //strcat(filename,".bin");
      //strcpy(fn1,filename);
      /*fp1 = fopen(fn01[j1], "rb");
      if (fp1==NULL)printf("File opening error\n");
      //strcpy(filename,"_jcnv");
      //strcat(filename,rankname);
      //strcat(filename,j1name);
      //strcat(filename,".bin");
      //strcpy(fn2,filename);
      fp2 = fopen(fn02[j1], "rb");
      if (fp2==NULL)printf("File opening error\n");
      //strcpy(filename,"_vav");
      //strcat(filename,rankname);
      //strcat(filename,j1name);
      //strcat(filename,".bin");
      //strcpy(fn3,filename);
      fp3 = fopen(fn03[j1], "rb");
      if (fp3==NULL)printf("File opening error\n");
      //ai= irnbs+j1irnbs;
      //aj= jcnbs+j1nz;
      //vals=valbs+j1nz;
      //nz=nzbs[j1];
      printf("rank %d j1nz %d nz %d\n",rank,j1nz,nz);
      //free(xi);
      //Multiply Bi by vi:
      if(VecSize-sumrowcolin<2000)l5=1;
      else l5=2000;
      int *rowi = (int*)calloc(l5*ha_ndblocks[j1+begblock[rank]],sizeof(int));
      int *coli = (int*)calloc(VecSize-sumrowcolin,sizeof(int));
      ha_cgetype *vali = (ha_cgetype*)calloc(l5*ha_ndblocks[j1+begblock[rank]],sizeof(ha_cgetype));
      fread(coli,sizeof(int),VecSize-sumrowcolin,fp2);
      fclose(fp2);
      //printf("rank %d irn %d jcn %d\n",rank,ha_ndblocks[j1+begblock[rank]],VecSize-sumrowcolin);
      if(VecSize-sumrowcolin<2000)for(i=0; i<VecSize-sumrowcolin; i++) {
          //if(rank==0)printf("rank %d i %d coli %d\n",rank,i,coli[i]);
          if(coli[i]<0)continue;
          fread(rowi,sizeof(int),coli[i],fp1);
          fread(vali,sizeof(ha_cgetype),coli[i],fp3);
          for(j=0; j<coli[i]; j++) {
            rowi[j]-=1;
            //printf("rank %d rowi %d\n",rank,rowi[j]);
            if(rowi[j]<ncolb-1) {
              for(j2=ai[rowi[j]]; j2<ai[rowi[j]+1]; j2++) {
                //printf("ib %d jb %d iv %d\n",j,aj[j2],iv);
                vecbivi[i+nrowc*aj[j2]]-=vali[j]*vals[j2];
              }
            } else {
              for(j2=ai[rowi[j]]; j2<nz; j2++) {
                vecbivi[i+nrowc*aj[j2]]-=vali[j]*vals[j2];
              }
            }
          }
        }
      else {
        l4=(VecSize-sumrowcolin)/l5;
        for(i=0; i<l4; i++) {
          //if(rank==0)printf("rank %d i %d coli %d\n",rank,i,coli[i]);
          j4=0;
          for(j3=0; j3<l5; j3++)if(coli[i*l5+j3]>=0)j4+=coli[i*l5+j3];
          if(j4==0)continue;
          fread(rowi,sizeof(int),j4,fp1);
          fread(vali,sizeof(ha_cgetype),j4,fp3);
          j4=0;
          for(j3=0; j3<l5; j3++) {
            j6=i*l5+j3;
            if(coli[j6]<0)continue;
            for(j=0; j<coli[j6]; j++) {
              rowi[j4]-=1;
              //printf("rank %d rowi %d\n",rank,rowi[j]);
              if(rowi[j4]<ncolb-1) {
                for(j2=ai[rowi[j4]]; j2<ai[rowi[j4]+1]; j2++) {
                  //printf("ib %d jb %d iv %d\n",j,aj[j2],iv);
                  vecbivi[j6+nrowc*aj[j2]]-=vali[j4]*vals[j2];
                }
              } else {
                for(j2=ai[rowi[j4]]; j2<nz; j2++) {
                  vecbivi[j6+nrowc*aj[j2]]-=vali[j4]*vals[j2];
                }
              }
              j4++;
            }
          }
        }
        j4=0;
        for(i=l4*l5; i<VecSize-sumrowcolin; i++)if(coli[i]>=0)j4+=coli[i];
        if(j4>0) {
          fread(rowi,sizeof(int),j4,fp1);
          fread(vali,sizeof(ha_cgetype),j4,fp3);
          j5=0;
          for(i=l4*l5; i<VecSize-sumrowcolin; i++) {
            //if(rank==0)printf("rank %d i %d coli %d\n",rank,i,coli[i]);
            if(coli[i]<0)continue;
            for(j=0; j<coli[i]; j++) {
              rowi[j5]-=1;
              //printf("rank %d rowi %d\n",rank,rowi[j]);
              if(rowi[j5]<ncolb-1) {
                for(j2=ai[rowi[j5]]; j2<ai[rowi[j5]+1]; j2++) {
                  //printf("ib %d jb %d iv %d\n",j,aj[j2],iv);
                  vecbivi[i+nrowc*aj[j2]]-=vali[j5]*vals[j2];
                }
              } else {
                for(j2=ai[rowi[j5]]; j2<nz; j2++) {
                  vecbivi[i+nrowc*aj[j2]]-=vali[j5]*vals[j2];
                }
              }
              j5++;
            }
          }
        }
      }
      printf("OKKK!\n");
      free(rowi);
      free(coli);
      free(vali);
      fclose(fp1);
      fclose(fp3);*/
      MatDestroy_SeqAIJ(submatBT);
      /*jcnv=realloc(jcnv,insize[8]*sizeof(forint));
      irnv=realloc(irnv,insize[8]*sizeof(forint));
      valuesv=realloc(valuesv,insize[8]*sizeof(ha_cgetype));
      printf("OKKKKK\n");
      MatCreateSeqAIJ(PETSC_COMM_SELF,ha_ndblocks[j1+rank*nmatinplus],VecSize-sumrowcolin,0,nnz,&submatv[j1]);
      //nnzmax=0;
      printf("rank %d nnz size %d col size %d\n",rank,ha_ndblocks[j1+rank*nmatinplus],VecSize-sumrowcolin);
      for(i=0; i<insize[8]-1; i++) {
        //if(valuesv[i]>1000000000)printf("i %d v %lf\n",i,valuesv[i]);
        //printf("r %d c %d v %lf nzz %d\n",irnv[i],jcnv[i],valuesv[i],nnz[irnv[i]]);
        //MatSetValues(submatv[j1],1,(PetscInt)&irnv[i],1,(PetscInt)&jcnv[i],&valuesv[i],INSERT_VALUES);
        MatSetValue(submatv[j1],irnv[i],jcnv[i],(PetscScalar)valuesv[i],INSERT_VALUES);
        //nnzmax+=nnzc[i];
      }
      printf("OK!!\n");
      MatAssemblyBegin(submatv[j1],MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(submatv[j1],MAT_FINAL_ASSEMBLY);*/
      //free(nnz);
      //free(irnv);
      //free(jcnv);
      //free(valuesv);
      //free(nnzc);
      printf("rank %d m %d n %d\n",rank,nrowb,ncolb);
      //VecGetSize(ui[j1],&nzv);
      //printf("vec n %d\n",nzv);
      //VecCreateSeq(PETSC_COMM_SELF,nrowb,&biui[j1]);
      //VecDuplicate(ui[j1],&biui[j1]);
      /*printf("nzv0 %d nzv %d rank %d OK????\n",nzv,insize[8],ha_ndblocks[j1+rank*nmatinplus]);
      //MatCreateSeqAIJ(PETSC_COMM_SELF,VecSize-sumrowcolin,ha_ndblocks[j1+rank*nmatinplus],PETSC_DEFAULT,NULL,&submatbivi[j1]);
      MatMatMult(submatB[j1],submatv[j1],MAT_INITIAL_MATRIX,PETSC_DEFAULT,&submatbivi);
      //if(rank==1){MatMatMult(submatA[j1],submatv[j1],MAT_INITIAL_MATRIX,PETSC_DEFAULT,&submatAA);
      //PetscViewerASCIIOpen(PETSC_COMM_SELF, "AsubAA.txt", &viewer);
      //ierr = MatView(submatAA,viewer);
      //CHKERRQ(ierr);
      //PetscViewerDestroy(&viewer);}
      aa=(Mat_SeqAIJ*)submatbivi->data;//*aa=subA->data;
      ai= aa->i;
      aj= aa->j;
      vals=aa->a;
      nz=aa->nz;
      nrow=submatbivi->rmap->n;
      ncol=submatbivi->cmap->n;
      printf("nz %d\n",nz);
      irn=(forint *) calloc (nz,sizeof(forint));
      if(nrow==1)for(i=0; i<nz; i++)irn[i]=0;
      else {
        if(ai[1]>0)for(j=ai[0]; j<ai[1]; j++)irn[j]=0;
        for(i=1; i<nrow-1; i++)if(ai[i+1]>ai[i]) {
            for(j=ai[i]; j<ai[i+1]; j++) {
              //printf("aici %d aicip1 %d j %d\n",aic[i],aic[i+1],j);
              irn[j]=i;
            }
          }
      }
      if(ai[nrow-1]<nz)for(j=ai[nrow-1]; j<nz; j++) {
          irn[j]=nrow-1;
        }
      /*for(i=0; i<nrow-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
          irn[j]=i;
        }
      for(j=ai[nrow-1]; j<nz; j++) {
        irn[j]=nrow-1;
      }*/
      //printf("i %d val %lf\n",i,vals[9452]);
      //printf("rank %d j1 %d vec %lf\n",rank,j1,vecbivi[22800]);
      //printf("0rank %d j1 %d vec %lf\n",rank,j1,vecbivi[16416]);
      /*for(i=0; i<nz; i++) {
        //if(rank==0&&j1==0)printf("rank %d j1 %d i %d aj %d irn %d indx %d vals %lf vecbivi %lf\n",rank,j1,i,aj[i],irn[i],aj[i]+ncol*irn[i],vals[i],vecbivi[aj[i]+ncol*irn[i]]);
        vecbivi[aj[i]+ncol*irn[i]]-=vals[i];
        //if(aj[i]+ncol*irn[i]==22800)printf("rank %d i %d val %lf\n",rank,i,vals[i]);
        //printf("rank %d j1 %d i %d aj %d irn %d indx %d vals %lf vecbivi %lf\n",rank,j1,i,aj[i],irn[i],aj[i]+ncol*irn[i],vals[i],vecbivi[aj[i]+ncol*irn[i]]);
      }*/
      //printf("1rank %d vec %lf\n",rank,vecbivi[16416]);
      //printf("1 rank %d j1 %d vec %lf\n",rank,j1,vecbivi[22800]);
      //free(irn);
      //MatDestroy_SeqAIJ(submatbivi);
      //MatMult(submatB[j1],ui1[j1],biui[j1]);
      //MatDestroy_SeqAIJ(submatB[j1]);
      //VecView(ui[j1],0);
      //VecGetArray(biui[j1],&vals);
      //for(i=0; i<nrowb; i++) {
      //VecGetValues(biui[j1],1,&i,&vecval);
      //vecbiui[i]+=vals[i];
      //}
      //VecRestoreArray(biui[j1],&vals);
      //VecDestroy(&biui[j1]);
      //free(insize);
      if(j1<0) {
        strcpy(tempchar,"zip -1 ");
        strcat(tempchar,fn01[j1]);
        strcat(tempchar,".zip ");
        strcat(tempchar,fn01[j1]);
        systemresult=system(tempchar);
        remove(fn01[j1]);
        strcpy(tempchar,"zip -1 ");
        strcat(tempchar,fn02[j1]);
        strcat(tempchar,".zip ");
        strcat(tempchar,fn02[j1]);
        systemresult=system(tempchar);
        remove(fn02[j1]);
        strcpy(tempchar,"zip -1 ");
        strcat(tempchar,fn03[j1]);
        strcat(tempchar,".zip ");
        strcat(tempchar,fn03[j1]);
        systemresult=system(tempchar);
        remove(fn03[j1]);
      }
      //j1nz+=nzbs[j1];
      //j1irnbs+=bns[j1];
      xi1indx+=ha_ndblocks[j1+begblock[rank]];
      time(&timeend);
      printf("Submatrix %d rank %d calculation time %f\n",j1,rank,difftime(timeend,timestr));
    }
    //MPI_Barrier(PETSC_COMM_WORLD);
  }
  //free(nzbs);//= realloc (nzbs,1*sizeof(uvadd));
  //free(bns);//= realloc (bns,1*sizeof(int));
  //free(bms);//= realloc (bms,1*sizeof(int));
  //free(irnbs);//= realloc (irnbs,1*sizeof(int));
  //free(jcnbs);//= realloc (jcnbs,1*sizeof(int));
  //free(valbs);//= realloc (valbs,1*sizeof(ha_cgetype));
  //ierr = PetscFree(biui);
  ierr = PetscFree(submatA);
  ierr = PetscFree(submatB);
  CHKERRQ(ierr);
  free(yi1);
  //MPI_Barrier(PETSC_COMM_WORLD);
  printf("Completed calculation of partitioned matrices! Rank %d\n",rank);
  timestr=clock();
  //forint vecbivisize=0;
  //if(rank==0) {
  vecbivisize=bivirowsize*bivicolsize;//(VecSize-sumrowcolin)*(VecSize-sumrowcolin);  //
  //printf("vecbivisize1 %d row %d col %d rank %d\n",vecbivisize,bivirowsize,bivicolsize,rank);
  //vecbiuisize=VecSize-sumrowcolin;
  //}
  //HaReduce(vecbivi,vecbivisize,mpisize,rank,mpisize-1);
  //HaReducebivi(vecbivi,bivinzrow,bivinzcol,bivinzrow0,bivinzcol0,bivinzrow1,bivinzcol1,bivirowsize,bivicolsize,mpisize,rank,mpisize-1);
  //int *biviindx1=(int *) calloc (1,sizeof(int));

  j=0;
  for(i=0; i<vecbivisize; i++)if(vecbivi[i]!=0)j++;
  unsigned int *biviindx1=(unsigned int *) calloc (j,sizeof(unsigned int));
  //ha_cgetype *vecbivi0= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
  int nz0,nz1,nz2,nz3=j;
  j2=0;
  for(i=0; i<vecbivisize; i++) {
    if(vecbivi[i]!=0) {
      vecbivi[j2]=vecbivi[i];
      j1=i/bivicolsize;
      j3=i-j1*bivicolsize;
      //printf("j1 %d j3 %d\n",i,bivinzcol[bivinzcol1[j3]]-1+(bivinzrow[bivinzrow1[j1]]-1)*bivicolsize);
      //biviindx[j2]=bivinzcol0[bivinzcol1[j3]]-1+(bivinzrow0[bivinzrow1[j1]]-1)*bivicolsize0;
      biviindx1[j2]=bivinzcol1[j3]+(unsigned int)bivinzrow1[j1]*nrowb;
      //if(biviindx1[j2]<0)printf("i %d j1 %d j3 %d biviindx1[j2] %d bivinzcol1[j3] %d bivinzrow1[j1] %d nrowb %d\n",i,j1,j3,biviindx1[j2],bivinzcol1[j3],bivinzrow1[j1],nrowb);
      //if(rank==0)printf("v %d\n",biviindx[j2]);
      j2++;
    }
  }
  vecbivi=realloc(vecbivi,j2*sizeof(ha_cgetype));
  nz0=j2;
  //printf("jjjj1 %d j2 %d rank %d\n",j,j2,rank);
  MPI_Reduce(&j,&j2,1, MPI_INT, MPI_MAX,mpisize-1,PETSC_COMM_WORLD);
  if(rank!=mpisize-1)j2=0;
  unsigned int *biviindx0= (unsigned int *) calloc (j2,sizeof(unsigned int));
  ha_cgetype *vecbivi0= (ha_cgetype *) calloc (j2,sizeof(ha_cgetype));
  //if(rank!=0)vecbivi=realloc(vecbivi,j*sizeof(ha_cgetype));

  for(j6=0; j6<mpisize-1; j6++) {
    j2=nz3;
    MPI_Bcast(&j2,1, MPI_INT,j6, PETSC_COMM_WORLD);
    j3=j6;
    MPI_Bcast(&j3,1, MPI_INT,j6, PETSC_COMM_WORLD);
    if(rank==j6) {
      if(SORD==1)MPI_Send(vecbivi,j2, MPI_DOUBLE,mpisize-1, 10, PETSC_COMM_WORLD);
      else MPI_Send(vecbivi,j2, MPI_FLOAT,mpisize-1, 10, PETSC_COMM_WORLD);
      MPI_Send(biviindx1,j2, MPI_INT,mpisize-1, 11, PETSC_COMM_WORLD);
    }
    //printf("j2 %d\n",j2);
    //printf("j1 %d rank000 %d\n",j6,rank);
    if(rank==mpisize-1) {
      //vecbivi0=realloc(vecbivi0,j2*sizeof(ha_cgetype));
      //biviindx=realloc(biviindx,j2*sizeof(int));
      //printf("j3 %d\n",j3);
      if(SORD==1)MPI_Recv(vecbivi0,j2, MPI_DOUBLE, j3, 10, PETSC_COMM_WORLD,&status);
      else MPI_Recv(vecbivi0,j2, MPI_FLOAT, j3, 10, PETSC_COMM_WORLD,&status);
      MPI_Recv(biviindx0,j2, MPI_INT, j3, 11, PETSC_COMM_WORLD,&status);
      //printf("in %d in0 %d j1 %d\n",biviindx1[0],biviindx0[0],j1);


      nz2=nz0;
      j=j2-1;
      j1=0;
      for(i=nz2-1; i>-1; i--) {
        if(biviindx1[i]>biviindx0[j]) {
          j1--;
          continue;
        }
        if(biviindx1[i]==biviindx0[j]) {
          j--;
          j1--;
          if(j==-1) {
            i--;
            break;
          }
          continue;
        }
        for(j3=j; j3>-1; j3--) {
          if(biviindx1[i]>biviindx0[j3]) {
            j1--;
            break;
          }
          if(biviindx1[i]==biviindx0[j3]) {
            j3--;
            j1--;
            break;
          }
          j1--;
        }
        j=j3;
        if(j==-1)break;
      }
      //printf("j1 %d j %d i %d v %d v0 %d\n",j1,j,i,biviindx1[i+1],biviindx0[j+1]);
      if(j>-1)j1-=j+1;//for(i=j;i>-1;i--){
      //j1--;
      //}
      //if(i==0)j1--;
      if(i>-1)j1-=i+1;//break
      //for(j=i;j>-1;j--){
      //j1--;
      //}
      nz1=-j1;
      /*j=0;
      nz1=0;
      nz2=nz0;
      printf("nz2 %d\n",nz2);
      j4=j2-1;
      printf("in %d in0 %d\n",biviindx1[0],biviindx0[0]);
      for(i=0; i<nz2; i++) {
        //printf("j %d i %d\n",j,i);
        if(biviindx1[i]<biviindx0[j]){nz1++;continue;}
        if(biviindx1[i]==biviindx0[j]){nz1++;j++;if(j==j4)break;continue;}
        for(j3=j;j3<j2;j3++){
          //printf("j3 %d\n",j3);
          if(biviindx1[i]<biviindx0[j3]){nz1++;break;}
          if(biviindx1[i]==biviindx0[j3]){nz1++;j3++;break;}
          nz1++;
        }
        j=j3;
        //nz1++;
        if(j==j4)break;
      }
      printf("i %d j %d nz1 %d\n",i,j,nz1);
      if(i<nz2-1)nz1+=nz2-1-i;
      if(j<j4)nz1+=j4-j;*/
      //printf("nz1 %d\n",nz1);
      vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
      //printf("hare nz1 %d\n",nz1);
      //if(vecbivi==NULL)printf("Relloc failed!");
      biviindx1=realloc(biviindx1,nz1*sizeof(unsigned int));
      //if(biviindx==NULL)printf("Relloc failed!");
      //printf("hare nz1 %d\n",nz1);
      j=j2-1;
      j1=nz1-1;
      for(i=nz2-1; i>-1; i--) {
        //printf("i %d j1 %d j %d v %d v0 %d v1 %d\n",i,j1,j,biviindx[i],biviindx0[j],biviindx[j1+1]);
        if(biviindx1[i]>biviindx0[j]) {
          vecbivi[j1]=vecbivi[i];
          biviindx1[j1]=biviindx1[i];
          j1--;
          continue;
        }
        if(biviindx1[i]==biviindx0[j]) {
          vecbivi[j1]=vecbivi[i]+vecbivi0[j];
          biviindx1[j1]=biviindx1[i];
          j--;
          j1--;
          if(j==-1)break;
          continue;
        }
        for(j3=j; j3>-1; j3--) {
          if(biviindx1[i]>biviindx0[j3]) {
            vecbivi[j1]=vecbivi[i];
            biviindx1[j1]=biviindx1[i];
            j1--;
            break;
          }
          if(biviindx1[i]==biviindx0[j3]) {
            vecbivi[j1]=vecbivi[i]+vecbivi0[j3];
            biviindx1[j1]=biviindx1[i];
            j3--;
            j1--;
            break;
          }
          biviindx1[j1]=biviindx0[j3];
          vecbivi[j1]=vecbivi0[j3];
          j1--;
        }
        j=j3;
        //biviindx[j1]=biviindx[i];
        //vecbivi[j1]=vecbivi[i];
        //j1--;
        if(j==-1)break;
      }
      //printf("j1 %d j %d i %d v %d v0 %d\n",j1,j,i,biviindx1[i+1],biviindx0[j+1]);
      if(j>-1)for(i=j; i>-1; i--) {
          biviindx1[j1]=biviindx0[i];
          vecbivi[j1]=vecbivi0[i];
          j1--;
        }
      //printf("j1 %d j %d i %d vm1 %d\n",j1,j,i,biviindx1[j1+1]);
      //for(i=0;i<2;i++)printf("vout %d\n",biviindx1[i]);
      nz0=nz1;


      //HaSpaVecAdd(vecbivi,biviindx1,nz,vecbivi0,biviindx0,j2);
      //printf("out in %d in0 %d j1 %d\n",biviindx1[0],biviindx0[0],j1);
    }
    //MPI_Barrier(PETSC_COMM_WORLD);
  }
  //printf("rank111 %d\n",rank);
  //free(biviindx1);
  free(biviindx0);
  free(vecbivi0);


  //HaReducebivi(vecbivi,nz0,biviindx1,bivinzrow,bivinzcol,bivinzrow0,bivinzcol0,bivinzrow1,bivinzcol1,bivirowsize,bivicolsize,nrowb,mpisize,rank,mpisize-1);
  //MatDestroyMatrices(nmatin,&submatB);
  //if(rank==mpisize-1)for(i=0; i<nz0; i++)printf("rank %d indx %d mat %lf\n",rank,biviindx[i],vecbivi[i]);
  //if(rank==0){
  /*forint vecbivisize=0,vecbiuisize=0;
  if(rank==0) {
    vecbivisize=(VecSize-sumrowcolin)*(VecSize-sumrowcolin);  //
    vecbiuisize=VecSize-sumrowcolin;
  }
  ha_cgetype *vecbivi0;//= (ha_cgetype *) calloc (vecbivisize,sizeof(ha_cgetype));
  //}
  //printf("rank %d vec %lf\n",rank,vecbivi[15055]);
  //MPI_Reduce(vecbivi,vecbivi0,(VecSize-sumrowcolin)*(VecSize-sumrowcolin), MPI_DOUBLE, MPI_SUM,0,PETSC_COMM_WORLD);
  //MPI_Reduce(vecbiui,vecbiui0,(VecSize-sumrowcolin), MPI_DOUBLE, MPI_SUM,0,PETSC_COMM_WORLD);
  //if(rank==0)MPI_Reduce(MPI_IN_PLACE,vecbivi,(VecSize-sumrowcolin)*(VecSize-sumrowcolin), MPI_DOUBLE, MPI_SUM,0,PETSC_COMM_WORLD);
  //else MPI_Reduce(vecbivi,vecbivi,(VecSize-sumrowcolin)*(VecSize-sumrowcolin), MPI_DOUBLE, MPI_SUM,0,PETSC_COMM_WORLD);
  MPI_Status   status;
  //MPI_Reduce(vecbivi,vecbivi0,(VecSize-sumrowcolin)*(VecSize-sumrowcolin), MPI_DOUBLE, MPI_SUM,0,PETSC_COMM_WORLD);
  j=0;
  if(rank!=0)for(i=0; i<(VecSize-sumrowcolin)*(VecSize-sumrowcolin); i++)if(vecbivi[i]!=0)j++;
  if(rank==0) j=0;//(VecSize-sumrowcolin)*(VecSize-sumrowcolin);
  //printf("jjjj %d rank %d\n",j,rank);
  int *biviindx= (int *) calloc (j,sizeof(int));
  j2=0;
  if(rank!=0)for(i=0; i<(VecSize-sumrowcolin)*(VecSize-sumrowcolin); i++) {
      if(vecbivi[i]!=0) {
        vecbivi[j2]=vecbivi[i];
        biviindx[j2]=i;
        j2++;
      }
    }
  //if(rank!=0)vecbivi=realloc(vecbivi,j*sizeof(ha_cgetype));
  //printf("jjjj1 %d rank %d\n",j,rank);

  for(j1=1; j1<mpisize; j1++) {
    j2=j;
    MPI_Bcast(&j2,sizeof(j2), MPI_BYTE,j1, PETSC_COMM_WORLD);
    j3=j1;
    MPI_Bcast(&j3,sizeof(j3), MPI_BYTE,j1, PETSC_COMM_WORLD);
    if(rank==j1) {
      MPI_Send(vecbivi,j, MPI_DOUBLE, 0, 0, PETSC_COMM_WORLD);
      MPI_Send(biviindx,j, MPI_INT, 0, 1, PETSC_COMM_WORLD);
    }
    //printf("j2 %d\n",j2);
    if(rank==0) {
      vecbivi0=realloc(vecbivi0,j2*sizeof(ha_cgetype));
      biviindx=realloc(biviindx,j2*sizeof(int));
      MPI_Recv(vecbivi0,j2, MPI_DOUBLE, j3, 0, PETSC_COMM_WORLD,&status);
      MPI_Recv(biviindx,j2, MPI_INT, j3, 1, PETSC_COMM_WORLD,&status);
      for(i=0; i<j2; i++) {
        //printf("biviindx %d\n",biviindx[i]);
        vecbivi[biviindx[i]]+=vecbivi0[i];
      }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  free(biviindx);
  free(vecbivi0);*/
  //free(bivinzrow0);
  //free(bivinzcol0);
  //printf("di ngu!\n");
  //if(rank==mpisize-1)printf("in1 %d\n",biviindx1[0]);
  if(rank!=mpisize-1) {
    //free(vecbivi);free(biviindx1);
    free(vecbivi);
    vecbivi=NULL;
    //vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
    free(biviindx1);
    biviindx1=NULL;
    //biviindx1=realloc(biviindx1,1*sizeof(unsigned int));
  }
  //printf("di ngu1!\n");
  //if(rank==mpisize-1)printf("in2 %d\n",biviindx1[0]);
  //if(rank==0)MPI_Reduce(MPI_IN_PLACE,vecbiui,(VecSize-sumrowcolin), MPI_DOUBLE, MPI_SUM,0,PETSC_COMM_WORLD);
  //else MPI_Reduce(vecbiui,vecbiui,(VecSize-sumrowcolin), MPI_DOUBLE, MPI_SUM,0,PETSC_COMM_WORLD);
  forint vecbiuisize=0;
  //if(rank==0) {
  //vecbivisize=(VecSize-sumrowcolin)*(VecSize-sumrowcolin);  //
  vecbiuisize=VecSize-sumrowcolin;
  //}
  //ha_cgetype *vecbiui0= (ha_cgetype *) calloc (vecbiuisize,sizeof(ha_cgetype));
  //MPI_Reduce(vecbiui,vecbiui0,(VecSize-sumrowcolin), MPI_DOUBLE, MPI_SUM,0,PETSC_COMM_WORLD);
  printf("after bivi rank %d\n",rank);
  HaReduce(vecbiui,vecbiuisize,mpisize,rank,mpisize-1);
  printf("Completed MPI_Reduce Operation! Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  //if(rank==0)printf("vec %lf\n",vecbivi0[15055]);
  //for(i=13528; i<13984; i++)printf("rank %d i %d vec %lf vecsum %lf\n",rank,i,vecbivi[i],vecbivi0[i]);
  if(rank!=mpisize-1)free(vecbiui);
  //if(rank!=0) {
  //free(vecbivi0);
  //free(vecbiui0);
  //free(vecbivi);
  //free(vecbiui);
  //}
  timestr=clock();
  ha_cgetype *xd;//= (ha_cgetype*)calloc(vecbiuisize,sizeof(ha_cgetype));
  ha_cgetype *x0;//,*vecbivi0;//= (ha_cgetype*)calloc(VecSize,sizeof(ha_cgetype));
  //x0=realloc (x0,1*sizeof(ha_cgetype));
  //int nz1,nz2;
  forint j7;
  if(rank==mpisize-1) {
    for(i=0; i<vecbiuisize; i++)yd[i]-=vecbiui[i];
    free(vecbiui);
    Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatD[0]->data;//*aa=subA->data;
    ai= aa->i;
    aj= aa->j;
    vals=aa->a;
    nz=aa->nz;
    nrow=submatD[0]->rmap->n;
    ncol=submatD[0]->cmap->n;
    //la=8;//5;
    printf("D ncol %d nrow %d nz %d\n",ncol,nrow,nz);
    /*forint *irn=(forint *) calloc (nz,sizeof(forint));
    if(nrow==1)for(i=0; i<nz; i++)irn[i]=0;
    else {
      if(ai[1]>0)for(j=ai[0]; j<ai[1]; j++)irn[j]=0;
      for(i=1; i<nrow-1; i++)if(ai[i+1]>ai[i]) {
          for(j=ai[i]; j<ai[i+1]; j++) {
            //printf("aici %d aicip1 %d j %d\n",aic[i],aic[i+1],j);
            irn[j]=i;
          }
        }
    }
    if(ai[nrow-1]<nz)for(j=ai[nrow-1]; j<nz; j++) {
        irn[j]=nrow-1;
      }*/
    /*for(i=0; i<nrow-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
        irn[j]=i;
      }
    for(j=ai[nrow-1]; j<nz; j++) {
      irn[j]=nrow-1;
    }*/
    //printf("vec %lf\n",vecbivi0[15055]);
    //ha_cgetype *valdleft=(ha_cgetype *) calloc (nz,sizeof(ha_cgetype));
    //forint *irndleft=(forint *) calloc (nz,sizeof(forint));
    //forint *jcndleft=(forint *) calloc (nz,sizeof(forint));
    unsigned int *biviindx0=(unsigned int *) calloc (nz,sizeof(unsigned int));
    j2=0;
    for(i=0; i<nrow-1; i++) {
      j7=(forint)i*ncol;
      for(j=ai[i]; j<ai[i+1]; j++) {
        biviindx0[j2]=aj[j]+j7;
        j2++;
      }
    }
    j7=(forint)i*ncol;
    for(j=ai[i]; j<nz; j++) {
      biviindx0[j2]=aj[j]+j7;
      //if(biviindx0[j2]<0)printf("biviindx0[j2] %d j %d j3 %d\n",biviindx0[j2],j,j3);
      j2++;
    }
    /*j2=0;
      for(i=0;i<nrow-1;i++){
        for(j=ai[i];j<ai[i+1];j++){
          if(bivinzrow[i]>0&&bivinzcol[aj[j]]>0)vecbivi[bivinzcol[aj[j]]+bivicolsize*(bivinzrow[i]-1)-1]+=vals[j];
          else{valdleft[j2]=vals[j];irndleft[j2]=i+1;jcndleft[j2]=aj[j]+1;j2++;}
        }
      }
        for(j=ai[i];j<nz;j++){
          if(bivinzrow[i]>0&&bivinzcol[aj[j]]>0)vecbivi[bivinzcol[aj[j]]+bivicolsize*(bivinzrow[i]-1)-1]+=vals[j];
          else{valdleft[j2]=vals[j];irndleft[j2]=i+1;jcndleft[j2]=aj[j]+1;j2++;}
        }*/

    j=nz-1;
    j1=0;
    printf("anz2 %d nz %d\n",nz2,nz);
    for(i=nz0-1; i>-1; i--) {
      //printf("i %d j1 %d j %d v %d v0 %d v1 %d\n",i,j1,j,biviindx[i],biviindx0[j],biviindx[j1+1]);
      if(biviindx1[i]>biviindx0[j]) {
        j1--;
        continue;
      }
      if(biviindx1[i]==biviindx0[j]) {
        j--;
        j1--;
        if(j==-1) {
          i--;
          break;
        }
        continue;
      }
      for(j3=j; j3>-1; j3--) {
        if(biviindx1[i]>biviindx0[j3]) {
          j1--;
          break;
        }
        if(biviindx1[i]==biviindx0[j3]) {
          j3--;
          j1--;
          break;
        }
        j1--;
      }
      j=j3;
      if(j==-1)break;
    }
    //printf("j1 %d j %d i %d v %d v0 %d\n",j1,j,i,biviindx1[i+1],biviindx0[j+1]);
    if(j>-1)j1-=j+1;//for(i=j;i>-1;i--){
    //j1--;
    //}
    //if(i==0)j1--;
    if(i>-1)j1-=i+1;//break
    //for(j=i;j>-1;j--){
    //j1--;
    //}
    nz1=-j1;

    /*j=0;
    nz1=0;
    nz2=nz0;
    printf("anz2 %d nz %d\n",nz2,nz);
    j4=nz-1;
    printf("in %d in0 %d\n",biviindx1[nz2-1],biviindx0[nz-1]);
    for(i=0; i<nz2; i++) {
      //printf("j %d i %d\n",j,i);
      if(biviindx1[i]<biviindx0[j]){nz1++;continue;}
      if(biviindx1[i]==biviindx0[j]){nz1++;j++;if(j==j4)break;continue;}
      for(j3=j;j3<nz;j3++){
        //printf("j3 %d\n",j3);
        if(biviindx1[i]<biviindx0[j3]){nz1++;break;}
        if(biviindx1[i]==biviindx0[j3]){nz1++;j3++;break;}
        nz1++;
      }
      j=j3;
      //nz1++;
      if(j==j4)break;
    }
    printf("i %d j %d nz1 %d\n",i,j,nz1);
    if(i<nz2-1)nz1+=nz2-1-i;
    if(j<j4)nz1+=j4-j;*/
    printf("nz1 %d\n",nz1);
    ldsize=ceil((laD/100.0)*nz1);
    vecbivi=realloc(vecbivi,ldsize*sizeof(ha_cgetype));
    //printf("hare nz1 %d\n",nz1);
    //if(vecbivi==NULL)printf("Relloc failed!");
    biviindx1=realloc(biviindx1,nz1*sizeof(unsigned int));
    //if(biviindx==NULL)printf("Relloc failed!");
    //printf("hare nz1 %d\n",nz1);
    j=nz-1;
    j1=nz1-1;
    printf("anz2 %d nz %d\n",nz0,nz);
    //for(i=0;i<10;i++)printf("vecbui20 %lf idx %d\n",vecbivi[i],biviindx1[i]);//
    //for(i=0;i<10;i++)printf("vecD20 %lf idx %d\n",vals[i],biviindx0[i]);//
    for(i=nz0-1; i>-1; i--) {
      //printf("i %d j1 %d j %d b %d\n",i,j1,j,biviindx1[j1-1]);
      if(biviindx1[i]>biviindx0[j]) {
        vecbivi[j1]=vecbivi[i];
        biviindx1[j1]=biviindx1[i];
        j1--;
        continue;
      }
      if(biviindx1[i]==biviindx0[j]) {
        vecbivi[j1]=vecbivi[i]+vals[j];
        biviindx1[j1]=biviindx1[i];
        j--;
        j1--;
        if(j==-1)break;
        continue;
      }
      for(j3=j; j3>-1; j3--) {
        if(biviindx1[i]>biviindx0[j3]) {
          vecbivi[j1]=vecbivi[i];
          biviindx1[j1]=biviindx1[i];
          j1--;
          break;
        }
        if(biviindx1[i]==biviindx0[j3]) {
          vecbivi[j1]=vecbivi[i]+vals[j3];
          biviindx1[j1]=biviindx1[i];
          j3--;
          j1--;
          break;
        }
        biviindx1[j1]=biviindx0[j3];
        vecbivi[j1]=vals[j3];
        j1--;
      }
      j=j3;
      //biviindx[j1]=biviindx[i];
      //vecbivi[j1]=vecbivi[i];
      //j1--;
      if(j==-1)break;
    }
    //for(i=0;i<nz1;i++)printf("i %d b %d v %lf\n",i,biviindx1[i],vecbivi[i]);
    //printf("j1 %d j %d i %d v %d v0 %d\n",j1,j,i,biviindx1[i+1],biviindx0[j+1]);
    if(j>-1)for(i=j; i>-1; i--) {
        biviindx1[j1]=biviindx0[i];
        vecbivi[j1]=vals[i];
        j1--;
      }
    //printf("j1 %d j %d i %d vm1 %d\n",j1,j,i,biviindx1[j1+1]);
    //for(i=0;i<2;i++)printf("vout %d\n",biviindx1[i]);
    nz0=nz1;


    //HaSpaVecAdd(vecbivi,biviindx,nz0,vals,biviindx0,nz);
    MatDestroy(&submatD[0]);
    free(biviindx0);
    //for(i=0;i<10;i++)printf("vecbui2 %lf idx %d\n",vecbivi[i],biviindx1[i]);//
    //for(i=0; i<nz; i++) {
    //printf("D irn %d jcn %d vals %lf i %d vec %lf\n",irn[i],aj[i],vals[i],aj[i]+ncol*irn[i],vecbivi0[aj[i]+ncol*irn[i]]);
    //vecbivi0[aj[i]+ncol*irn[i]]+=vals[i];
    //vecbivi[aj[i]+ncol*irn[i]]+=vals[i];
    //printf("D irn %d jcn %d vals %lf vec %lf\n",irn[i],aj[i],vals[i],vecbivi0[aj[i]+ncol*irn[i]]);
    //}
    //printf("vec1 %lf\n",vecbivi0[15055]);
    //free(irn);
    nz=nz0;
    //for(i=0; i<vecbivisize; i++)if(vecbivi[i]!=0)nz++;
    //nz+=j2;
    //nz=VecSize-sumrowcolin;
    //printf("vecbivi ncol %d nrow %d nz %d\n",VecSize-sumrowcolin,ncol,nz);
    //irn=(forint *) calloc (nz,sizeof(forint));
    int *irn1=(int *) calloc (ldsize,sizeof(int));
    int *jcn=(int *) calloc (ldsize,sizeof(int));
    //ha_cgetype *values= (ha_cgetype *) calloc (la*nz,sizeof(ha_cgetype));
    for(i=0; i<nz; i++) {
      j3=biviindx1[i]/ncol;
      irn1[i]=j3+1;
      jcn[i]=biviindx1[i]-j3*ncol+1;
      //if(irn1[i]<0)printf("i %d bivi %d\n",i,biviindx1[i]);
      //printf("irn %d jcn %d biviindx %d j3 %d col %d\n",irn1[i],jcn[i],biviindx1[i],j3,ncol);
      //values[i]=vecbivi[i];
      //printf("j1 %d vec %lf\n",j1,vecbivi[i]);
    }
    /*for(i=0; i<vecbivisize; i++)if(vecbivi[i]!=0) {
        j3=i/bivicolsize;
        irn1[j]=bivinzrow1[j3]+1;
        jcn[j]=bivinzcol1[i-j3*bivicolsize]+1;
        values[j]=vecbivi[i];
        //printf("irn %d jcn %d values %lf\n",irn1[j],jcn[j],values[j]);
        j++;
      }
    for(i=0; i<j2; i++){
        irn1[j]=irndleft[i];
        jcn[j]=jcndleft[i];
        values[j]=valdleft[i];
        //printf("irn1 %d jcn1 %d values1 %lf\n",irn1[j],jcn[j],values[j]);
        j++;
      }*/
    //free(valdleft);
    //free(irndleft);
    //free(jcndleft);
    free(bivinzrow);
    free(bivinzcol);
    free(bivinzrow1);
    free(bivinzcol1);
    //free(biviindx1);
    //for(i=0;i<j;i++)printf("jcn %d\n",jcn[i]);
    //for(i=0; i<nz; i++){
    //irn[i]=i+1;
    //jcn[i]=i+1;
    //values[i]=1;
    //}
    //free(vecbivi);
    //vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
    free(biviindx1);
    biviindx1=NULL;
    //biviindx1=realloc(biviindx1,1*sizeof(unsigned int));
    int *insizeD=(int *) calloc (4,sizeof(int));
    insizeD[0]=nrow;
    insizeD[1]=ncol;
    insizeD[2]=nz;
    //if(1>=(VecSize-sumrowcolin)*(VecSize-sumrowcolin)/nz)
    insizeD[3]=laD;
    printf("ok\n");
    xd=(ha_cgetype *) calloc (vecbiuisize,sizeof(ha_cgetype));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
    printf("ok!!!!!!!!\n");
    spec48_ssol2la_(insizeD,irn1,jcn,vecbivi,yd,xd);
    free(yd);
    x0=(ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));//realloc (x0,*sizeof(ha_cgetype));
    for(i=0; i<vecbiuisize; i++) {
      x0[indicesC[i]]=xd[i];
      //printf("x1d %lf\n",xd[i]);
    }
    free(irn1);
    free(jcn);
    //free(values);
    free(insizeD);
  }
  free(biviindx1);
  free(vecbivi);
  free(indicesC);
  if(rank!=mpisize-1) {
    free(bivinzrow);
    free(bivinzcol);
    free(bivinzrow1);
    free(bivinzcol1);
    xd=(ha_cgetype *) calloc (vecbiuisize,sizeof(ha_cgetype));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
    x0=(ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));//realloc (x0,*sizeof(ha_cgetype));
  }
  PetscFree(submatD);
  //if(rank!=0) {
  //free(vecbiui0); //free(vecbivi0);
  //}
  MPI_Barrier(PETSC_COMM_WORLD);
  if(SORD==1)MPI_Bcast(xd, vecbiuisize, MPI_DOUBLE,mpisize-1, PETSC_COMM_WORLD);
  else MPI_Bcast(xd, vecbiuisize, MPI_FLOAT,mpisize-1, PETSC_COMM_WORLD);
  printf("Interface problem solution Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  //VecCreateSeq(PETSC_COMM_SELF,VecSize-sumrowcolin,&vecxd);
  //PetscInt *idexix=(PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
  //for(i=0; i<(VecSize-sumrowcolin); i++)idexix[i]=i;
  //VecSetValues(vecxd,VecSize-sumrowcolin,idexix,xd, INSERT_VALUES);
  //if(rank==1)VecView(vecxd,0);
  //free(xd);
  //free(idexix);
  timestr=clock();
  xi1indx=0;
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      ha_cgetype *biui0= (ha_cgetype *) calloc (ha_ndblocks[j1+begblock[rank]],sizeof(ha_cgetype));
      //if(j1<10)strcpy(j1name,"000");
      //if(j1<100&&j1>9)strcpy(j1name,"00");
      //if(j1<1000&&j1>99)strcpy(j1name,"0");
      //sprintf(filename, "%d",j1);
      //strcat(j1name,filename);
      //strcpy(filename,"_irnv");
      //strcat(filename,rankname);
      //strcat(filename,j1name);
      //strcat(filename,".bin");
      //strcpy(fn1,filename);
      //strcpy(filename,"_jcnv");
      //strcat(filename,rankname);
      //strcat(filename,j1name);
      //strcat(filename,".bin");
      //strcpy(fn2,filename);
      //strcpy(filename,"_vav");
      //strcat(filename,rankname);
      //strcat(filename,j1name);
      //strcat(filename,".bin");
      //strcpy(fn3,filename);

      if(j1<0) {
        strcpy(tempchar,"unzip ");
        strcat(tempchar,fn01[j1]);
        strcat(tempchar,".zip");
        systemresult=system(tempchar);
        strcpy(tempchar,fn01[j1]);
        strcat(tempchar,".zip");
        remove(tempchar);
        strcpy(tempchar,"unzip ");
        strcat(tempchar,fn02[j1]);
        strcat(tempchar,".zip");
        systemresult=system(tempchar);
        strcpy(tempchar,fn02[j1]);
        strcat(tempchar,".zip");
        remove(tempchar);
        strcpy(tempchar,"unzip ");
        strcat(tempchar,fn03[j1]);
        strcat(tempchar,".zip");
        systemresult=system(tempchar);
        strcpy(tempchar,fn03[j1]);
        strcat(tempchar,".zip");
        remove(tempchar);
      }

      fp1 = fopen(fn01[j1], "rb");
      printf("fn1 %s\n",fn01[j1]);
      if (fp1==NULL)printf("File opening error\n");
      fp2 = fopen(fn02[j1], "rb");
      if (fp2==NULL)printf("File opening error\n");
      fp3 = fopen(fn03[j1], "rb");
      if (fp3==NULL)printf("File opening error\n");
      la1=ceil((insize[j1*insizes+9]/100.0)*insize[j1*insizes+2]);
      int *irne = (int*)calloc(la1,sizeof(int));
      int *keep = (int*)calloc(insize[j1*insizes+12],sizeof(int));
      ha_cgetype *vale = (ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
      freadresult=fread(irne,sizeof(int),la1,fp1);
      freadresult=fread(keep,sizeof(int),insize[j1*insizes+12],fp2);
      freadresult=fread(vale,sizeof(ha_cgetype),la1,fp3);
      fclose(fp1);
      fclose(fp2);
      fclose(fp3);
      remove(fn01[j1]);
      remove(fn02[j1]);
      remove(fn03[j1]);
      //for(i=0;i<la;i++)printf("vav %f\n",vale[i]);
      Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatC[j1]->data;//*aa=subA->data;
      //printf("0 rank %d j1 %d vec %lf\n",rank,j1,vecbivi[22800]);
      ai= ac->i;
      aj= ac->j;
      vals=ac->a;
      nz=ac->nz;
      nrow=submatC[j1]->rmap->n;
      ncol=submatC[j1]->cmap->n;
      printf("nrow %d\n",nrow);
      ha_cgetype *be0 = (ha_cgetype*)calloc(nrow,sizeof(ha_cgetype));
      for(i=0; i<nrow-1; i++) {
        for(j=ai[i]; j<ai[i+1]; j++) {
          be0[i]+=vals[j]*xd[aj[j]];
        }
      }
      for(j=ai[i]; j<nz; j++) {
        be0[i]+=vals[j]*xd[aj[j]];
      }
      MatDestroy(&submatC[j1]);
      spec48_esol_(insize+j1*insizes,irne,vale,keep,be0,biui0);
      printf("here!\n");
      free(irne);
      free(keep);
      free(vale);
      free(be0);
      ISGetIndices(colindices[j1],&nindices);
      xi1point=xi1+xi1indx;
      for(i=0; i<ha_ndblocks[j1+begblock[rank]]; i++) {
        x0[nindices[i]]+=xi1point[i]-biui0[i];//vals[i]-biui0[i];//vecval;
      }
      ISRestoreIndices(colindices[j1],&nindices);
      free(biui0);
      xi1indx+=ha_ndblocks[j1+begblock[rank]];
    }
    //MPI_Barrier(PETSC_COMM_WORLD);
  }
  free(xi1);
  free(insize);
  ierr = PetscFree(submatC);
  //free(vnz);
  free(begblock);
  printf("Solution calculation Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  //printf("Ok0 rank %d\n",rank);
  MPI_Barrier(PETSC_COMM_WORLD);
  //printf("Ok0 rank %d\n",rank);
  free(xd);
  timestr=clock();
  if(SORD==1)MPI_Allreduce(x0,x1,VecSize, MPI_DOUBLE, MPI_SUM,PETSC_COMM_WORLD);
  else MPI_Allreduce(x0,x1,VecSize, MPI_FLOAT, MPI_SUM,PETSC_COMM_WORLD);
  printf("Reduce solution Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  //MPI_Bcast(x1, VecSize*sizeof(ha_cgetype), MPI_BYTE,0, PETSC_COMM_WORLD);
  //VecDestroy(&vecxd);
  //printf("Ok1 rank %d\n",rank);
  free(x0);
  free(fn01);
  free(fn02);
  free(fn03);
  //printf("Ok2 rank %d\n",rank);
  //for(i=0; i<VecSize; i++)if(rank==0)printf("x1 %lf\n",x1[i]);
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&colindices[i]);
    CHKERRQ(ierr);
    //VecDestroy(&ui1[i]);
  }
  PetscFree(colindices);
  //printf("Ok rank %d\n",rank);
  //PetscFree(ui1);
  printf("OK rui rank %d\n",rank);
  return 0;
}

int HaDBBDParSol(Mat A, Vec b, ha_cgetype *x1, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, uvadd ndblock,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laD,PetscReal cntl3) {//,bool iter
  IS *rowindices,*colindices,*Cindices,*Bindices,*BBindices,*BBiindices;
  const PetscInt *nindices;
  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nrowcolin,sumrowcolin,i,i1,j,j0,j1,j2,j3,j4,j5,j6,l0,l1,l2,l3,l4,l5,rank1,proc1=0,nnzmax,j1nz,j1irnbs;
  Mat *submatA,*submatC,*submatB,*submatD,submatBT,submatCT,*submatBB;
  Vec vecxd;
  PetscInt *ai,*aj,*aic,*ajc;
  PetscInt nrow,ncol,nz,nrowc,ncolc,nrowb,ncolb,nzc,nzv;
  uvadd lasize,ldsize;
  //forint iv,jv,la1;
  MPI_Status   status;
  clock_t timestr,timeend;
  //forint forintwrp;
  int j2int,la1;
  size_t freadresult;
  //long long int nrow,ncol,nz,nrowc,ncolc,nzc,nzv;
  PetscScalar *vals,*valsc,vecval;
  ha_cgetype vval;
  PetscErrorCode ierr;
  PetscViewer viewer;
  FILE* fp1,*fp2,*fp3;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  char filename[1024],rankname[1024],j1name[1024],tempchar[1024];//,fn1[1024],fn2[1024],fn3[1024]
  if(rank<10)strcpy(rankname,"000");
  if(rank<100&&rank>9)strcpy(rankname,"00");
  if(rank<1000&&rank>99)strcpy(rankname,"0");
  if(rank>=1000)rankname[0]='\0';
  sprintf(filename, "%d",rank);
  strcat(rankname,filename);
  nmatin=(uvadd)ndblock/mpisize;
  int nmatinBB=1;
  if(nmatin*mpisize<ndblock)nmatinplus=nmatin+1;
  else nmatinplus=nmatin;
  for(i=0; i<mpisize; i++)if(rank+1<=ndblock-mpisize*nmatin)nmatin++;
  //forint *vnz= (forint *) calloc (nmatin,sizeof(forint));
  //printf("rank %d nmatin %d\n",rank,nmatin);
  begblock[rank]=nmatin;
  for(i=0; i<mpisize; i++) {
    j=rank;
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    j=j2;
    MPI_Bcast(&begblock[j2],1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  j2=begblock[0];
  begblock[0]=0;
  for(i=1; i<mpisize; i++) {
    j1=begblock[i];
    begblock[i]=begblock[i-1]+j2;
    j2=j1;
  }
  //j2=begblock[0];
  //for(i=0; i<mpisize; i++)begblock[i]-=j2;
  //if(rank==1)for(i=0; i<mpisize; i++)printf("begblock %d\n",begblock[i]);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&rowindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&colindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&Cindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&Bindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatinBB*sizeof(IS **),&BBindices);
  CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatinBB*sizeof(IS **),&BBiindices);
  //CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatA);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatC);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatB);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatinBB*sizeof(Mat *),&submatBB);
  CHKERRQ(ierr);
  PetscScalar **yi1= (PetscScalar**)calloc(nmatin,sizeof(PetscScalar*));
  for (i=0; i<nmatin; i++) yi1[i] = (PetscScalar*)calloc(ha_ndblocks[i+begblock[rank]],sizeof(PetscScalar));
  sumrowcolin=0;
  for(j=0; j<ndblock; j++)sumrowcolin+=ha_ndblocks[j];
  PetscInt *indicesC= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
  PetscInt *indicesB= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
  //PetscScalar *yd = (PetscScalar*)calloc(VecSize-sumrowcolin,sizeof(PetscScalar));
  PetscInt *offblock= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
  PetscInt *offblockrow= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
  ha_cgetype *vecbiui= (ha_cgetype *) calloc (VecSize-sumrowcolin,sizeof(ha_cgetype));
  offblock[0]=0;
  for(j=0; j<ndblock; j++)offblock[j+1]=offblock[j]+countvarintra1[j+1]-countvarintra1[j]-ha_ndblocks[j];
  offblockrow[0]=0;
  for(j=0; j<ndblock; j++)offblockrow[j+1]=offblockrow[j]+counteqnoadd[j]-ha_ndblocks[j];
  //printf("OK123\n");

  char **fn01= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn01[i] = (char*)calloc(1024,sizeof(char));
  char **fn02= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn02[i] = (char*)calloc(1024,sizeof(char));
  char **fn03= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn03[i] = (char*)calloc(1024,sizeof(char));

  timestr=clock();
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      if(j1<10)strcpy(j1name,"000");
      if(j1<100&&j1>9)strcpy(j1name,"00");
      if(j1<1000&&j1>99)strcpy(j1name,"0");
      if(j1>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j1);
      strcat(j1name,filename);
      strcpy(filename,"_irnv");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn01[j1],filename);
      strcpy(filename,"_keep");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn02[j1],filename);
      strcpy(filename,"_vav");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn03[j1],filename);
    }
  }
  printf("OK012345 rank %d mpis %d\n",rank,mpisize);
  for(j=0; j<mpisize; j++) {
    printf("j %d\n",j);
    for(j3=0; j3<nmatinplus; j3++) {
      if(j3==nmatin)i=j3-1;
      else i=j3;
      //printf("rank %d OK1!!! i %d j %d nmatin %d\n",rank,i,j,nmatin);
      bfirst=counteq[i+begblock[rank]];
      bend=ha_ndblocks[i+begblock[rank]];
      MPI_Bcast(&bfirst,1, MPI_INT,j, PETSC_COMM_WORLD);
      MPI_Bcast(&bend,1, MPI_INT,j, PETSC_COMM_WORLD);
      PetscInt *indices= (PetscInt *) calloc (bend,sizeof(PetscInt));
      ha_cgetype *yi0 = (ha_cgetype*)calloc(bend,sizeof(ha_cgetype));
      for(j1=0; j1<bend; j1++) { //ha_ndblocks[i+begblock[j]]
        j2=ha_rows[bfirst+j1]+bfirst;
        if(j2>=Istart&&j2<Iend)indices[j1]=j2;
        else indices[j1]=-1;
      }
      //printf("rank %d OK1!!! i %d j %d\n",rank,i,j);
      VecGetValues(b,bend,indices,yi0);
      //VecView(b,0);
      //for(i=0;i<bend;i++)printf("rank %d yi0 %lf\n",rank,yi0[i]);
      //printf("rank %d OK2!!! i %d j %d nmatin %d\n",rank,i,j,nmatin);

      HaReduce(yi0,(forint)bend,mpisize,rank,j);
      if(rank==j&&i<nmatin) {
        memcpy (yi1[i],yi0,bend*sizeof(ha_cgetype));
      }
      free(indices);
      free(yi0);
    }
  }
  printf("OK12345 rank %d\n",rank);
  j1=0;
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (ha_ndblocks[i+begblock[rank]],sizeof(PetscInt));
    //printf("ha_row 0 %d 588 %d\n",ha_rows[0],ha_rows[588]);
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_rows[bfirst+j]+bfirst;
      //if(rank==0&&i==0)printf("indices %d rank %d\n",indices[j],ha_ndblocks[i+rank*nmatinplus]);
    }

    bend=counteq[i+begblock[rank]]+counteqnoadd[i+begblock[rank]];
    for(j=offblockrow[i+begblock[rank]]; j<offblockrow[i+1+begblock[rank]]; j++) {
      j1=j-offblockrow[i+begblock[rank]];
      indicesB[j]=ha_rows[bend-j1-1]+bfirst;
      //printf("rank %d rindices %d bend %d j1 %d j %d\n",rank,indicesB[j],ha_rows[bend+j1],j1,j);
    }
    //printf("i %d\n",i);
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
    //printf("i %d\n",i);
    //ISView(rowindices[i],PETSC_VIEWER_STDOUT_WORLD);
    bfirst=countvarintra1[i+begblock[rank]];
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_cols[bfirst+j]+bfirst;
    }
    bend=countvarintra1[i+1+begblock[rank]];
    for(j=offblock[i+begblock[rank]]; j<offblock[i+1+begblock[rank]]; j++) {
      j1=j-offblock[i+begblock[rank]];
      indicesC[j]=ha_cols[bend-j1-1]+bfirst;
      //printf("rank %d cindices %d bend %d j1 %d j %d\n",rank,indicesC[j],ha_cols[bend-j1-1],j1,j);
    }

    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,colindices+i);
    free(indices);
    //ISView(colindices[i],PETSC_VIEWER_STDOUT_WORLD);
    //if(ha_ndblocks[i+rank*nmatinplus]!=0)yi0=realloc (yi0,1*sizeof(PetscScalar));
  }
  printf("0okkkk\n");
  j1=offblock[ndblock];
  for(i=countvarintra1[ndblock]; i<VecSize; i++) {
    indicesC[j1]=i;
    //printf("rank %d indices %d j1 %d\n",rank,indicesC[j1],j1);
    j1++;
  }
  printf("1okkkk\n");
  j1=offblockrow[ndblock];
  for(i=counteq[ndblock-1]+counteqnoadd[ndblock-1]; i<VecSize; i++) {
    indicesB[j1]=i;
    //if(rank==0)printf("rank %d indices %d j1 %d\n",rank,indicesB[j1],j1);
    j1++;
  }
  printf("2okkkk\n");
  //for(i=0; i<VecSize-sumrowcolin; i++)if(rank==0)printf("0rank %d bindx %d\n",rank,indicesB[i]);
  for(i=0; i<mpisize; i++) {
    j=offblockrow[begblock[rank]];
    j1=offblockrow[nmatin+begblock[rank]]-offblockrow[begblock[rank]];
    //printf("rank %d i %d j %d j1 %d\n",rank,i,j,j1);
    MPI_Bcast(&j,1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&j1,1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&indicesB[j],j1, MPI_INT,i, PETSC_COMM_WORLD);
    //if(rank==i)VecGetValues(b,j1,&indicesB[j],&yd[j]);
  }
  //printf("rank0 %d\n",rank);
  //for(i=0; i<offblockrow[ndblock]; i++)if(rank==1)printf("1rank %d bindx %d\n",rank,indicesB[i]);
  //for(i=0; i<VecSize-sumrowcolin; i++)if(rank==2)printf("indiceB %d\n",indicesB[i]);
  printf("rank1 %d\n",rank);
  //for(i=0; i<VecSize-sumrowcolin; i++){if(rank==1)printf("i %d\n",i);VecGetValues(b,1,&indicesB[i],&yd[i]);}
  //for(i=0; i<mpisize; i++)if(rank==i)

  for(i=0; i<VecSize-sumrowcolin; i++)if(indicesB[i]>=Istart) {
      j=i;
      break;
    }
  for(j2=i; j2<VecSize-sumrowcolin; j2++)if(indicesB[j2]>=Iend) {
      j1=j2;
      break;
    }
  if(j2==VecSize-sumrowcolin)j1=VecSize-sumrowcolin;
  printf("rank %d j %d j1 %d istart %d iend %d\n",rank,j,j1,Istart,Iend);
  VecGetValues(b,j1-j,&indicesB[j],&vecbiui[j]);//implicite yd save mem
  ierr = VecDestroy(&b);
  CHKERRQ(ierr);
  //MPI_Barrier(PETSC_COMM_WORLD);
  //if(rank==mpisize-1)VecGetValues(b,VecSize-sumrowcolin-bend,indicesB+bend,yd+bend);
  printf("rank2 %d\n",rank);
  //MPI_Reduce(yd,yd0,VecSize-sumrowcolin, MPI_DOUBLE, MPI_SUM,0,PETSC_COMM_WORLD);
  //forintwrp=VecSize-sumrowcolin;
  //HaReduce(vecbiui,forintwrp,mpisize,rank,mpisize-1);
  //if(rank!=mpisize-1)free(yd);
  //if(rank!=0)free(yd0);
  //printf("Vec %d sumrow %d\n",VecSize,sumrowcolin);
  //for(i=0; i<VecSize-sumrowcolin; i++)if(rank==0)printf("rank %d yd %f\n",rank,yd0[i]);
  for(i=0; i<mpisize; i++) {
    j=offblock[begblock[rank]];
    j1=offblock[nmatin+begblock[rank]]-offblock[begblock[rank]];
    //printf("rank %d i %d j %d j1 %d\n",rank,i,j,j1);
    MPI_Bcast(&j,1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&j1,1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&indicesC[j],j1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  printf("OK132\n");
  for(i=0; i<nmatin; i++) {
    ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesC,PETSC_COPY_VALUES,Cindices+i);
    ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesB,PETSC_COPY_VALUES,Bindices+i);
  }
  free(indicesB);
  //for(i=0; i<VecSize-sumrowcolin; i++)printf("rank %d indiceB %d\n",rank,indicesB[i]);
  //for(i=0; i<VecSize-sumrowcolin; i++)printf("rank %d indiceC %d\n",rank,indicesC[i]);
  free(offblock);
  free(offblockrow);
  PetscInt *indexBB= (PetscInt *) calloc (VecSize,sizeof(PetscInt));
  for(i=0; i<VecSize; i++)indexBB[i]=i;
  ISCreateGeneral(PETSC_COMM_SELF,VecSize,indexBB,PETSC_COPY_VALUES,BBindices);
  printf("Begin partitioning rank %d\n",rank);
  //if(rank==0)ISView(Bindices[0],PETSC_VIEWER_STDOUT_SELF);
  //if(rank==0)ISView(Cindices[0],PETSC_VIEWER_STDOUT_SELF);
  //if(rank==0)ISView(rowindices[0],PETSC_VIEWER_STDOUT_SELF);
  //if(rank==0)ISView(colindices[0],PETSC_VIEWER_STDOUT_SELF);
  ierr = MatCreateSubMatrices(A,nmatin,rowindices,colindices,MAT_INITIAL_MATRIX,&submatA);
  printf("A\n");
  ierr = MatCreateSubMatrices(A,nmatin,rowindices,Cindices,MAT_INITIAL_MATRIX,&submatC);
  printf("C\n");
  //ierr = MatCreateSubMatrices(A,nmatin,Bindices,colindices,MAT_INITIAL_MATRIX,&submatB);
  ierr = MatCreateSubMatrices(A,nmatinBB,Bindices,BBindices,MAT_INITIAL_MATRIX,&submatBB);
  ierr = MatDestroy(&A);
  CHKERRQ(ierr);
  //MatView(submatBB[0],0);
  PetscInt BBrow,BBcol,CCrow,CCcol;
  MatGetSize(submatBB[0],&BBrow,&BBcol);
  MatGetSize(submatC[0],&CCrow,&CCcol);
  printf("bbrow %d bbcol %d ccrow %d cccol %d\n",BBrow,BBcol,BBcol,CCcol);
  PetscInt *indexBBi= (PetscInt *) calloc (BBrow,sizeof(PetscInt));
  for(i=0; i<BBrow; i++)indexBBi[i]=i;
  for(i=0; i<nmatin; i++) {
    ISCreateGeneral(PETSC_COMM_SELF,BBrow,indexBBi,PETSC_COPY_VALUES,Bindices+i);
  }
  //ISCreateGeneral(PETSC_COMM_SELF,BBrow,indexBBi,PETSC_COPY_VALUES,BBiindices);
  ierr = MatCreateSubMatrices(submatBB[0],nmatin,Bindices,colindices,MAT_INITIAL_MATRIX,&submatB);
  printf("B\n");
  if(rank==mpisize-1)proc1=1;
  ierr = PetscMalloc(proc1*sizeof(Mat *),&submatD);
  CHKERRQ(ierr);
  printf("D\n");
  //MatCreateSubMatrices(A,proc1,Bindices,Cindices,MAT_INITIAL_MATRIX,&submatD);
  MatCreateSubMatrices(submatBB[0],proc1,Bindices,Cindices,MAT_INITIAL_MATRIX,&submatD);
  printf("End partitionning rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  //CHKERRQ(ierr);
  //PetscViewerASCIIOpen(PETSC_COMM_SELF, "A.txt", &viewer);
  //PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_DENSE );
  //ierr = MatView(A,viewer);
  //CHKERRQ(ierr);
  //PetscViewerDestroy(&viewer);
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&rowindices[i]);
    CHKERRQ(ierr);
    //ierr = ISDestroy(colindices+i);
    //CHKERRQ(ierr);
    ierr = ISDestroy(&Cindices[i]);
    CHKERRQ(ierr);
    ierr = ISDestroy(&Bindices[i]);
    CHKERRQ(ierr);
  }
  ierr = ISDestroy(&BBindices[0]);
  //ierr = ISDestroy(&BBiindices[0]);
  PetscFree(rowindices);
  PetscFree(Cindices);
  PetscFree(Bindices);
  PetscFree(BBindices);
  //PetscFree(BBiindices);
  PetscFree(indexBB);
  PetscFree(indexBBi);
  ierr = MatDestroy(&submatBB[0]);
  PetscFree(submatBB);
  CHKERRQ(ierr);
  //PetscViewerASCIIOpen(PETSC_COMM_SELF, "b.txt", &viewer);
  //ierr = VecView(b,0);
  //CHKERRQ(ierr);
  //PetscViewerDestroy(&viewer);
//  if(rank<0) {
//    //PetscViewerASCIIOpen(PETSC_COMM_SELF, "Asub.txt", &viewer);
//    //PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_DENSE );
//    //ierr = MatView(submatA[0],viewer);
//    //CHKERRQ(ierr);
//    //PetscViewerDestroy(&viewer);
//    PetscViewerASCIIOpen(PETSC_COMM_SELF, "C.txt", &viewer);
//    //PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_DENSE );
//    ierr = MatView(submatC[0],viewer);
//    //CHKERRQ(ierr);
//    PetscViewerDestroy(&viewer);
//    PetscViewerASCIIOpen(PETSC_COMM_SELF, "B.txt", &viewer);
//    //PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_DENSE );
//    ierr = MatView(submatB[0],viewer);
//    //CHKERRQ(ierr);
//    PetscViewerDestroy(&viewer);
//    PetscViewerASCIIOpen(PETSC_COMM_SELF, "D.txt", &viewer);
//    //PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_DENSE );
//    ierr = MatView(submatD[0],viewer);
//    //CHKERRQ(ierr);
//    PetscViewerDestroy(&viewer);
//  }
  printf("OK1 rank %d sumrowcolin %d\n",rank,sumrowcolin);
  ha_cgetype *xi1 = (ha_cgetype*)calloc(sumrowcolin,sizeof(ha_cgetype));
  int insizes=17;
  int *insize=(int *) calloc (insizes*nmatin,sizeof(int));
  int *bivinzrow=(int *) calloc (VecSize-sumrowcolin,sizeof(int));
  int *bivinzcol=(int *) calloc (VecSize-sumrowcolin,sizeof(int));
  PetscReal *ccolnorms=(PetscReal *) calloc (BBrow,sizeof(PetscReal));
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      //MatTranspose(submatC[j1],MAT_INITIAL_MATRIX,&submatBT);
      //Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatBT->data;//*aa=subA->data;
      //ai= ac->i;
      //nz=ac->nz;
      //ncol=submatC[j1]->cmap->n;
      //for(i=0; i<ncol-1; i++)if(ai[i]!=ai[i+1])bivinzcol[i]++;
      //if(ai[i]<nz)bivinzcol[i]++;
      //MatDestroy(&submatBT);
      MatGetColumnNorms(submatC[j1],NORM_1,ccolnorms);
      for(i=0; i<BBrow; i++)if(ccolnorms[i]>0)bivinzcol[i]++;
      Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatB[j1]->data;//*aa=subA->data;
      ai= ab->i;
      nz=ab->nz;
      nrow=submatB[j1]->rmap->n;
      for(i=0; i<nrow-1; i++)if(ai[i]!=ai[i+1])bivinzrow[i]++;
      if(ai[i]<nz)bivinzrow[i]++;
    }
  }
  free(ccolnorms);
  int bivirowsize=1,bivicolsize=1;
  long int vecbivisize=0,li,lj;
  for(i=0; i<VecSize-sumrowcolin; i++) {
    if(bivinzrow[i]>0) {
      bivinzrow[i]=bivirowsize;
      bivirowsize++;
    }
    if(bivinzcol[i]>0) {
      bivinzcol[i]=bivicolsize;
      bivicolsize++;
    }
  }
  bivirowsize--;
  bivicolsize--;
  long int *bivinzrow1=(long int *) calloc (bivirowsize,sizeof(long int));
  long int *bivinzcol1=(long int *) calloc (bivicolsize,sizeof(long int));
  j=0;
  j2=0;
  for(li=0; li<VecSize-sumrowcolin; li++) {
    if(bivinzrow[li]>0) {
      bivinzrow1[j]=li*BBrow;
      //if(bivinzrow1[j]<0)printf("WARNINGS!!!! Index out of range!\n");
      j++;
    }
    if(bivinzcol[li]>0) {
      bivinzcol1[j2]=li;
      j2++;
    }
    bivinzrow[li]=(bivinzrow[li]-1)*bivicolsize;//faster MSOL
  }
  vecbivisize=((long int)bivirowsize)*((long int)bivicolsize);
  if(vecbivisize<0)printf("Overflow detected, the interface problem is too big!!!!!!\n");
  ha_cgetype *vecbivi= (ha_cgetype *) calloc (vecbivisize,sizeof(ha_cgetype));
  printf("vecbivisize %ld rank %d\n",vecbivisize,rank);

  ha_cgetype *xi1point;
  uvadd xi1indx=0;
  int jthrd,nthrd=1;
  //clock_t timestr1;
  nthrd=omp_get_max_threads();
  uvadd *nthrds= (uvadd *) calloc (nthrd+1,sizeof(uvadd));
  uvadd *nthrds1= (uvadd *) calloc (nthrd+1,sizeof(uvadd));
  #pragma omp parallel private(jthrd,j1,j,timestr,submatCT,ai,aj,vals,nz,nrow,ncol,lasize,aic,ajc,valsc,nzc,nrowc,ncolc,nrowb,ncolb,xi1point,xi1indx,i,timeend) shared(submatA,insize,nthrds,nthrds1)
  {
  jthrd=omp_get_thread_num();
  int windx=0,bindx,eindx;
  windx=90/nthrd;
  bindx=windx*jthrd;
  eindx=bindx+windx;
  windx=bindx;
  #pragma omp for
  for(j1=0; j1<nmatinplus; j1++){
    if(j1<nmatin) {
      //#pragma omp atomic
      nthrds[jthrd+1]+=ha_ndblocks[j1+begblock[rank]];
    }
    //#pragma omp atomic
    nthrds1[jthrd+1]+=1;
  }
  //printf("OK1 nj0 %ld rank %ld jth %d\n",nthrds1[jthrd],nthrds1[jthrd+1],jthrd);
  if(jthrd==0){
  for(j1=1; j1<nthrd+1; j1++){
    nthrds[j1]+=nthrds[j1-1];
    nthrds1[j1]+=nthrds1[j1-1];
    //printf("OK1 j1 %d rank %d jth %d nj1 %ld\n",j1,rank,jthrd,nthrds1[j1]);
  }
  }
  #pragma omp barrier
  xi1indx=nthrds[jthrd];
  //printf("OK1 nj1 %ld rank %ld jth %d\n",nthrds1[jthrd],nthrds1[jthrd+1],jthrd);
  //printf("jth %d OK1 j1 %ld j1end %ld xi %ld\n",jthrd,nthrds1[jthrd],nthrds1[jthrd+1],xi1indx);
  //#pragma omp for
//       #pragma omp critical
//       {
  for(j1=nthrds1[jthrd]; j1<nthrds1[jthrd+1]; j1++) {
    //printf("OK1 j1 %d rank %d jth %d\n",j1,rank,jthrd);
    if(j1<nmatin) {
      time(&timestr);//=clock();
      //printf("OK1 j1 %d rank %d\n",j1,rank);
      Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatA[j1]->data;//*aa=subA->data;
      //printf("OK1 j1 %d\n",j1);
      MatTranspose(submatC[j1],MAT_INITIAL_MATRIX,&submatCT);
      //MatDestroy_SeqAIJ(submatC[j1]);
      Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatCT->data;//*aa=subA->data;
      //printf("0 rank %d j1 %d vec %lf\n",rank,j1,vecbivi[0]);
      ai= aa->i;
      aj= aa->j;
      vals=aa->a;
      nz=aa->nz;
      nrow=submatA[j1]->rmap->n;
      ncol=submatA[j1]->cmap->n;
      //printf("nz %d\n",nz);
      //la=2;
      lasize=ceil((laA/100.0)*nz);
      int *irn=(int *) calloc (lasize,sizeof(int));
      int *jcn=(int *) calloc (lasize,sizeof(int));
      ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
      //ha_cgetype *xi = (ha_cgetype*)calloc(ha_ndblocks[j1+begblock[rank]],sizeof(ha_cgetype));
      //int *insize=(int *) calloc (12,sizeof(int));
      //long long int *insize=(long long int *) calloc (10,sizeof(long long int));
      insize[j1*insizes+13]=bivirowsize;
      insize[j1*insizes+14]=bivicolsize;
      insize[j1*insizes]=nrow;
      insize[j1*insizes+1]=ncol;
      insize[j1*insizes+2]=nz;
      aic= ac->i;
      ajc= ac->j;
      valsc=ac->a;
      nzc=ac->nz;
      //for(i=0;i<nzc;i++)++ajc[i];
      nrowc=submatCT->rmap->n;
      ncolc=submatCT->cmap->n;
      //printf("nrowc %d ncolc %d nzc %d\n",nrowc,ncolc,nzc);
      //int *irnc=(int *) calloc (nrowc,sizeof(int));
      //int *irnc1=(int *) calloc (nrowc,sizeof(int));
      //int *jcnc=(int *) calloc (nzc,sizeof(int));
      //ha_cgetype *valuesc= (ha_cgetype *) calloc (nzc,sizeof(ha_cgetype));
      insize[j1*insizes+3]=nrowc;
      insize[j1*insizes+4]=ncolc;
      //for(i=0; i<nrow; i++)if(rank==1)printf("0aici %d\n",ai[i]);
      /*if(nrow==1)for(i=0; i<nz; i++)irn[i]=1;//jcn[i]=aj[i]+1;}
      else {
        if(ai[1]>0)for(j=ai[0]; j<ai[1]; j++)irn[j]=1;//jcn[j]=aj[j]+1;}
        for(i=1; i<nrow-1; i++)if(ai[i+1]>ai[i]) {
            for(j=ai[i]; j<ai[i+1]; j++) {
              irn[j]=i+1;
              //jcn[j]=aj[j]+1;
              //if(rank==1&&(i==5||i==6))printf("aici %d aicip1 %d j %d irn %d\n",aic[i],aic[i+1],j,irn[j]);
            }
          }
      }
      if(ai[nrow-1]<nz)for(j=ai[nrow-1]; j<nz; j++) {
          irn[j]=nrow;
          //jcn[j]=aj[j]+1;
        }*/
      //for(i=0; i<nrow-1; i++)for(j=ai[i]; j<ai[i+1]; j++)irn[j]=i+1;
      //for(j=ai[i]; j<nz; j++)irn[j]=nrow;
      //j=0;
      //for(i=0; i<nz; i++) {
      //if(vals[i]!=0) {
      //irn[i]=irn[i];
      //jcn[i]=aj[i]+1;
      //values[i]=vals[i];
      //j++;
      //}
      //printf("irn %d jcn %d val %lf\n",irn[i],jcn[i],vals[i]);
      //}
      memcpy (irn+nz,ai,(nrow+1)*sizeof(PetscInt));
      memcpy (jcn,aj,nz*sizeof(PetscInt));
      memcpy (values,vals,nz*sizeof(ha_cgetype));
      //printf("nz %d j %d\n",nz,j);
      insize[j1*insizes+2]=nz;
      //for(i=0; i<nrowc-1; i++)if(rank==1)printf("0aici %d\n",irnc[i]);
      //for(i=0; i<nrowc-1; i++)irnc[i]=aic[i+1]-aic[i];//printf("i %d irnc %d\n",i,irnc[i]);}
      //irnc[i]=nzc-aic[i];//printf("i %d irnc %d\n",i,irnc[i]);
      //for(i=0; i<nzc; i++) {
      //  jcnc[i]=ajc[i]+1;
      //printf("jnc %d valsc %lf\n",jcnc[i],valsc[i]);
      //}
      //memcpy (valuesc,valsc,nzc*sizeof(ha_cgetype));
      insize[j1*insizes+5]=nzc;
      MatDestroy(&submatA[j1]);
      insize[j1*insizes+9]=laA;
      insize[j1*insizes+16]=lasize;
      insize[j1*insizes+10]=rank;
      insize[j1*insizes+11]=j1;
      //MatGetSize(submatB[j1],&nrowb,&ncolb);
      nrowb=submatB[j1]->rmap->n;
      ncolb=submatB[j1]->cmap->n;
      //printf("rank1 %d nzv size %d sizemat %d\n",rank,nzv,sizeof(valsc));
      //MatTranspose(submatB[j1],MAT_INITIAL_MATRIX,&submatBT);
      //MatDestroy(&submatB[j1]);
      Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatB[j1]->data;//*aa=subA->data;
      ai= ab->i;//irnbs+j1irnbs;
      aj= ab->j;//jcnbs+j1nz;
      vals=ab->a;//valbs+j1nz;
      nz=ab->nz;//nzbs[j1];
      insize[j1*insizes+6]=nrowb;//ncolb;
      insize[j1*insizes+7]=ncolb;//nrowb;
      insize[j1*insizes+8]=nz;
      insize[j1*insizes+15]=7+windx;//7+rank+jthrd+windx;//(j4+jthrd)%90+7;
      windx++;
      if(windx==eindx)windx=bindx;
      //for(i=0;i<nz;i++)++aj[i];
      xi1point=xi1+xi1indx;
      //for(i=0;i<xi1indx;i++)xi1[i]=0;
      //for(i=0;i<30524;i++)xi1point[i]=0;
      //for(i=0;i<10;i++)printf("vecbui0 %lf\n",vecbivi[i]);//bivirowsize*bivicolsize
      spec48m_msol_(insize+j1*insizes,irn,jcn,values,yi1[j1],xi1point,aic,ajc,valsc,ai,aj,vals,vecbivi,bivinzrow,bivinzcol);
      //printf("OK1 j1 %d rank %d\n",j1,rank);
      MatDestroy(&submatCT);
      //for(i=0;i<10;i++)printf("vecbui1 %lf\n",vecbivi[i]);//
      //printf("rank %d nzv1 %d OK???? xi1 %d\n",rank,nzv,xi1point-xi1);
      //if(rank==0)for(i=0;i<ncolb;i++)printf("rank %d i %d xi %lf\n",rank,i,xi1point[i]);
      //vnz[j1]=insize[8];
      free(yi1[j1]);
      //free(irnc);
      //free(jcnc);
      //free(valuesc);
      free(irn);
      free(jcn);
      free(values);
      //Multiply Bi by ui:
      for(i=0; i<nrowb-1; i++) {
        for(j=ai[i]; j<ai[i+1]; j++) {
          #pragma omp atomic
          vecbiui[i]-=vals[j]*xi1point[aj[j]];//xi[i];
        }
      }
      for(j=ai[i]; j<nz; j++) {
        #pragma omp atomic
        vecbiui[i]-=vals[j]*xi1point[aj[j]];//xi[i];
      }
      MatDestroy(&submatB[j1]);//submatBT);
      //printf("rank %d m %d n %d\n",rank,nrowb,ncolb);
//      if(j1<0) {
//        strcpy(tempchar,"zip -1 ");
//        strcat(tempchar,fn01[j1]);
//        strcat(tempchar,".zip ");
//        strcat(tempchar,fn01[j1]);
//        system(tempchar);
//        remove(fn01[j1]);
//        strcpy(tempchar,"zip -1 ");
//        strcat(tempchar,fn02[j1]);
//        strcat(tempchar,".zip ");
//        strcat(tempchar,fn02[j1]);
//        system(tempchar);
//        remove(fn02[j1]);
//        strcpy(tempchar,"zip -1 ");
//        strcat(tempchar,fn03[j1]);
//        strcat(tempchar,".zip ");
//        strcat(tempchar,fn03[j1]);
//        system(tempchar);
//        remove(fn03[j1]);
//      }
      //j1nz+=nzbs[j1];
      //j1irnbs+=bns[j1];
      xi1indx+=ha_ndblocks[j1+begblock[rank]];
      time(&timeend);
      printf("Submatrix %d rank %d thrd %d calculation time %f\n",j1,rank,jthrd,difftime(timeend,timestr));
    }
//       }
    //MPI_Barrier(PETSC_COMM_WORLD);
  }
  }
  free(nthrds);
  free(nthrds1);
  free(bivinzcol);
  free(bivinzrow);
  ierr = PetscFree(submatA);
  ierr = PetscFree(submatB);
  CHKERRQ(ierr);
  free(yi1);
  //MPI_Barrier(PETSC_COMM_WORLD);
  printf("Completed calculation of partitioned matrices! Rank %d\n",rank);
  timestr=clock();
  //vecbivisize=bivirowsize*bivicolsize;//(VecSize-sumrowcolin)*(VecSize-sumrowcolin);  //
  //halfvec=vecbivisize/2;

  #pragma omp parallel private(li) reduction(+:lj)
  {
  lj=0;
  #pragma omp for
  //for(li=vecbivisize-1; li!=-1; --li)if(vecbivi[li]!=0)lj++;
  for(li=0; li<vecbivisize; ++li){
    if(vecbivi[li]!=0)lj++;
  }
  }
  long int *obiviindx1,*biviindx1,*obiviindx0,*biviindx0;
  //obiviindx1=realloc(obiviindx1,sizeof(long int));
  //obiviindx0=realloc(obiviindx0,sizeof(long int));
  //obiviindx0[0]=-1;
  obiviindx1=realloc(obiviindx1,(lj+1)*sizeof(long int));
  obiviindx1[0]=-1;
  //int *obiviindx1=(int *) calloc (j,sizeof(int));
  //obiviindx1[0]=-1;
  biviindx1=obiviindx1;
  biviindx1++;
  //ha_cgetype *vecbivi0= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
  long int nz0,nz1,nz2,nz3=lj,halfj2;
  //my_vec_comz_(vecbivi,biviindx1,bivinzcol1,bivinzrow1,&bivicolsize,&vecbivisize,&j2int);
  //j2=j2int;
  long int lj2=0;
  for(li=0; li<vecbivisize; ++li) {
    if(vecbivi[li]!=0) {
      vecbivi[lj2]=vecbivi[li];
      //j1=i/bivicolsize;
      //j3=i-j1*bivicolsize;
      //printf("j1 %d j3 %d\n",i,bivinzcol[bivinzcol1[j3]]-1+(bivinzrow[bivinzrow1[j1]]-1)*bivicolsize);
      //biviindx[j2]=bivinzcol0[bivinzcol1[j3]]-1+(bivinzrow0[bivinzrow1[j1]]-1)*bivicolsize0;
      biviindx1[lj2]=bivinzcol1[li%bivicolsize]+bivinzrow1[li/bivicolsize];//[i-j1*bivicolsize]+bivinzrow1[j1];//*nrowb;
      //if(biviindx1[j2]<halfvec)halfj2=j2;
      //if(biviindx1[j2]<0)printf("i %ld j1 %d j3 %d biviindx1[j2] %ld bivinzcol1[j3] %ld bivinzrow1[j1] %ld nrowb %d\n",li,j1,j3,biviindx1[j2],bivinzcol1[j3],bivinzrow1[j1],nrowb);
      //if(rank==0)printf("v %d\n",biviindx[j2]);
      ++lj2;
    }
  }
  free(bivinzrow1);
  free(bivinzcol1);
  vecbivi=realloc(vecbivi,lj2*sizeof(ha_cgetype));
  nz0=lj2;
  printf("jjjj1 %ld j2 %ld rank %d nz0 %ld nz3 %ld\n",lj,lj2,rank,nz0,nz3);
  //MPI_Reduce(&j,&j2,1, MPI_INT, MPI_MAX,mpisize-1,PETSC_COMM_WORLD);
  //if(rank!=mpisize-1)j2=0;
  //unsigned int *biviindx0= (unsigned int *) calloc (j2,sizeof(unsigned int));
  //ha_cgetype *vecbivi0= (ha_cgetype *) calloc (j2,sizeof(ha_cgetype));
  //if(rank!=0)vecbivi=realloc(vecbivi,j*sizeof(ha_cgetype));
  /*unsigned int *biviindx0= (unsigned int *) calloc (1,sizeof(unsigned int));
  ha_cgetype *vecbivi0= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
  for(j6=0; j6<mpisize-1; j6++) {
    j2=halfj2;//nz3;
    nz2=nz0-halfj2;
    MPI_Bcast(&j2,1, MPI_INT,j6, PETSC_COMM_WORLD);
    j3=j6;
    MPI_Bcast(&j3,1, MPI_INT,j6, PETSC_COMM_WORLD);
    j4=nz2;//nz3;
    MPI_Bcast(&j4,1, MPI_INT,mpisize-1, PETSC_COMM_WORLD);
    if(rank==j6) {
      MPI_Send(vecbivi,j2, MPI_DOUBLE,mpisize-1, 10, PETSC_COMM_WORLD);
      MPI_Send(biviindx1,j2, MPI_INT,mpisize-1, 11, PETSC_COMM_WORLD);
      memmove (vecbivi,vecbivi+halfj2,nz2*sizeof(ha_cgetype));
      memmove (biviindx1,biviindx1+halfj2,nz2*sizeof(unsigned int));
      biviindx1=realloc(biviindx1,nz2*sizeof(unsigned int));
      vecbivi=realloc(vecbivi,nz2*sizeof(ha_cgetype));
      biviindx0=realloc(biviindx0,j4*sizeof(unsigned int));
      vecbivi0=realloc(vecbivi0,j4*sizeof(ha_cgetype));
      //unsigned int *biviindx0= (unsigned int *) calloc (j4,sizeof(unsigned int));
      //ha_cgetype *vecbivi0= (ha_cgetype *) calloc (j4,sizeof(ha_cgetype));
      MPI_Recv(vecbivi0,j4, MPI_DOUBLE,mpisize-1, 14, PETSC_COMM_WORLD,&status);
      MPI_Recv(biviindx0,j4, MPI_INT,mpisize-1, 15, PETSC_COMM_WORLD,&status);
      if(biviindx1[nz2-1]<biviindx0[j4-1]) {
        j=0;
        //j1=0;
        i1=0;
        for(i=0; i<nz2; i++) {
          while (biviindx1[i]>biviindx0[j])j++;//i1++;}
          //if(biviindx1[i]==biviindx0[j])j1++;
          if(biviindx1[i]<biviindx0[j])i1++;
        }
        nz1=j4+i1;//nz2+j2-j1;
        //printf("nz1 %d nz2 %d j1 %d j2 %d j %d i1 %d\n",nz1,nz2,j1,j2,j,i1);
      } else {
        j=0;
        //j1=0;
        i1=0;
        for(i=0; i<j4; i++) {
          while (biviindx0[i]>biviindx1[j]) j++;//i1++;}
          //if(i<10)printf("i %d j %d x0 %d x1 %d",i,j,biviindx0[i],biviindx1[j]);
          if(biviindx0[i]<biviindx1[j])i1++;
          //if(biviindx0[i]==biviindx1[j])j1++;
        }
        nz1=nz2+i1;//nz2+j2-j1;
      }
      vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
      biviindx1=realloc(biviindx1,nz1*sizeof(unsigned int));
      j=j4-1;
      j1=nz2-1;
      i=nz1-1;
      while(i>-1) {
        if(biviindx1[j1]==biviindx0[j]) {
            vecbivi[i]=vecbivi[j1]+vecbivi0[j];
            biviindx1[i]=biviindx1[j1];
            j1--;
            j--;
            i--;
            continue;
          }
        i1=0;
        while(biviindx1[j1]>biviindx0[j]) {
          vecbivi[i]=vecbivi[j1];
          biviindx1[i]=biviindx1[j1];
          i1++;
          j1--;
          i--;
        }
        while(biviindx1[j1]<biviindx0[j]) {
            vecbivi[i]=vecbivi0[j];
            biviindx1[i]=biviindx0[j];
            j--;
            i--;
          }
      }
      printf("!!!!!!!!!!!!rank %d i %d j %d j1 %d\n",rank,i,j,j1);
      biviindx0=realloc(biviindx1,sizeof(unsigned int));
      vecbivi0=realloc(vecbivi,sizeof(ha_cgetype));
      j5=nz1;
      //MPI_Bcast(&j5,1, MPI_INT,j6, PETSC_COMM_WORLD);
      MPI_Send(&j5,1, MPI_INT,mpisize-1, 16, PETSC_COMM_WORLD);
      printf("!!!!!!!!!!!!i %d j %d j1 %d j5 %d\n",i,j,j1,j5);
      MPI_Send(vecbivi,j5, MPI_DOUBLE,mpisize-1, 17, PETSC_COMM_WORLD);
      MPI_Send(biviindx1,j5, MPI_INT,mpisize-1, 18, PETSC_COMM_WORLD);
      biviindx1=realloc(biviindx1,sizeof(unsigned int));
      vecbivi=realloc(vecbivi,sizeof(ha_cgetype));
      //nz0=nz1;
    }
    if(rank==mpisize-1) {
      biviindx0=realloc(biviindx0,j2*sizeof(unsigned int));
      vecbivi0=realloc(vecbivi0,j2*sizeof(ha_cgetype));
      //unsigned int *biviindx0= (unsigned int *) calloc (j2,sizeof(unsigned int));
      //ha_cgetype *vecbivi0= (ha_cgetype *) calloc (j2,sizeof(ha_cgetype));
      MPI_Recv(vecbivi0,j2, MPI_DOUBLE, j3, 10, PETSC_COMM_WORLD,&status);
      MPI_Recv(biviindx0,j2, MPI_INT, j3, 11, PETSC_COMM_WORLD,&status);
      MPI_Send(vecbivi+halfj2,j4, MPI_DOUBLE,j3, 14, PETSC_COMM_WORLD);
      MPI_Send(biviindx1+halfj2,j4, MPI_INT,j3, 15, PETSC_COMM_WORLD);
      biviindx1=realloc(biviindx1,halfj2*sizeof(unsigned int));
      vecbivi=realloc(vecbivi,halfj2*sizeof(ha_cgetype));
      //nz2=nz0;
      if(biviindx1[halfj2-1]<biviindx0[j2-1]) {
        j=0;
        //j1=0;
        i1=0;
        for(i=0; i<halfj2; i++) {
          while (biviindx1[i]>biviindx0[j])j++;//i1++;}
          //if(biviindx1[i]==biviindx0[j])j1++;
          if(biviindx1[i]<biviindx0[j])i1++;
        }
        nz1=j2+i1;//nz2+j2-j1;
        //printf("nz1 %d nz2 %d j1 %d j2 %d j %d i1 %d\n",nz1,nz2,j1,j2,j,i1);
      } else {
        j=0;
        //j1=0;
        i1=0;
        for(i=0; i<j2; i++) {
          while (biviindx0[i]>biviindx1[j]) j++;//i1++;}
          //if(i<10)printf("i %d j %d x0 %d x1 %d",i,j,biviindx0[i],biviindx1[j]);
          if(biviindx0[i]<biviindx1[j])i1++;
          //if(biviindx0[i]==biviindx1[j])j1++;
        }
        nz1=halfj2+i1;//nz2+j2-j1;
      }
      vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
      biviindx1=realloc(biviindx1,nz1*sizeof(unsigned int));
      //printf("OKKK!\n");
      j=j2-1;
      j1=halfj2-1;
      i=nz1-1;
      while(i>-1) {
        if(biviindx1[j1]==biviindx0[j]) {
            vecbivi[i]=vecbivi[j1]+vecbivi0[j];
            biviindx1[i]=biviindx1[j1];
            j1--;
            j--;
            i--;
            continue;
          }
        i1=0;
        while(biviindx1[j1]>biviindx0[j]) {
          vecbivi[i]=vecbivi[j1];
          biviindx1[i]=biviindx1[j1];
          i1++;
          j1--;
          i--;
        }
        while(biviindx1[j1]<biviindx0[j]) {
            vecbivi[i]=vecbivi0[j];
            biviindx1[i]=biviindx0[j];
            j--;
            i--;
          }
      }
      printf("!!!!!!!!!!!!i %d j %d j1 %d\n",i,j,j1);
      halfj2=nz1;
      MPI_Recv(&j5,1, MPI_INT, j3, 16, PETSC_COMM_WORLD,&status);
      printf("OKKKK! j5 %d nz1 %d\n",j5,nz1);
      biviindx0=realloc(biviindx0,j5*sizeof(unsigned int));
      vecbivi0=realloc(vecbivi0,j5*sizeof(ha_cgetype));
      printf("OKKKK!\n");
      //unsigned int *biviindx0= (unsigned int *) calloc (j5,sizeof(unsigned int));
      //ha_cgetype *vecbivi0= (ha_cgetype *) calloc (j5,sizeof(ha_cgetype));
      MPI_Recv(vecbivi0,j5, MPI_DOUBLE, j3, 17, PETSC_COMM_WORLD,&status);
      printf("OKKKK!\n");
      MPI_Recv(biviindx0,j5, MPI_INT, j3, 18, PETSC_COMM_WORLD,&status);
      nz0=nz1+j5;
      biviindx1=realloc(biviindx1,nz0*sizeof(unsigned int));
      vecbivi=realloc(vecbivi,nz0*sizeof(ha_cgetype));
      memcpy (vecbivi+nz1, vecbivi0, j5*sizeof(ha_cgetype) );
      vecbivi0=realloc(vecbivi0,sizeof(ha_cgetype));
      memcpy (biviindx1+nz1, biviindx0, j5*sizeof(unsigned int) );
      biviindx0=realloc(biviindx0,sizeof(unsigned int));
    }
  }
  free(vecbivi0);
  free(biviindx0);*/
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len,name_len_max,name_beg,class_size,color,group_size,ha_id,group_size1,group_size11,ha_id1;
//   MPI_Comm HA_COMM,HA1_COMM;
//   MPI_Get_processor_name(processor_name, &name_len);
//   //printf("rank %d name len %d proc name %s\n",rank,name_len,processor_name);
//   MPI_Allreduce(&name_len,&name_len_max,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
//   name_len_max++;
//   char *vec_pr_name=(char *) calloc (mpisize*name_len_max,sizeof(char));
//   char *vec_pr_sname=(char *) calloc (mpisize*name_len_max,sizeof(char));
//   name_beg=rank*name_len_max;
//   for(i=name_beg; i<name_len+name_beg; i++)vec_pr_name[i]=processor_name[i-name_beg];
//   vec_pr_name[i]='\0';
//   //printf("rank %d name len %d proc name %s\n",rank,name_len,vec_pr_name+name_beg);
//   MPI_Barrier(PETSC_COMM_WORLD);
//   MPI_Allreduce(vec_pr_name,vec_pr_sname,mpisize*name_len_max,MPI_CHAR,MPI_SUM,PETSC_COMM_WORLD);
//   //if(rank==0)for(i=0;i<mpisize;i++)printf("snamessss %s\n",vec_pr_sname+i*name_len_max);
//   for(i=0; i<name_len_max; i++)vec_pr_name[i]=vec_pr_sname[i];
//   j=1;
//   for(i=0; i<mpisize; i++) {
//     for(j1=0; j1<j; j1++) {
//       if(strcmp(vec_pr_sname+i*name_len_max,vec_pr_name+j1*name_len_max)==0) {
//         break;
//       }
//     }
//     if(j1==j) {
//       //printf("rank %d sname %s name %s\n",rank,vec_pr_sname+i*name_len_max,vec_pr_name+j1*name_len_max);
//       for(j6=0; j6<name_len_max; j6++)vec_pr_name[j*name_len_max+j6]=vec_pr_sname[j1*name_len_max+j6];
//       j++;
//     }
//   }
//   class_size=j;
//   //printf("class size %d\n",class_size);
//   //if(rank==0)for(i=0;i<class_size;i++)printf("sname %s\n",vec_pr_sname+i*name_len_max);
//   for(i=0; i<class_size; i++) {
//     if(strcmp(processor_name,vec_pr_name+i*name_len_max)==0) {
//       color=i;
//       break;
//     }
//   }
//   free(vec_pr_name);
//   free(vec_pr_sname);
//   MPI_Comm_split(PETSC_COMM_WORLD,color,rank,&HA_COMM);
  MPI_Comm_rank( HA_COMM, &ha_id);
  MPI_Comm_size(HA_COMM,&group_size);
  //printf("rank %d color %d group %d ha_id %d\n",rank,color,group_size,ha_id);
  if(ha_id==group_size-1)color=1;
  else color=0;
//   MPI_Comm_split(PETSC_COMM_WORLD,color,rank,&HA1_COMM);
  MPI_Comm_rank( HA1_COMM, &ha_id1);
  MPI_Comm_size(HA1_COMM,&group_size1);
  group_size11=group_size1;
  MPI_Bcast(&group_size11,1, MPI_INT,mpisize-1, PETSC_COMM_WORLD);
  //printf("rank %d color %d group11 %d\n",rank,color,group_size11);
  int sindx01,sindx02,sindx03;

  if(ha_id!=group_size-1) {
    lj2=nz3;
    //printf("rank %d color %d group11 %d\n",rank,color,group_size11);
    MPI_Send(&lj2,1, MPI_LONG,group_size-1, 10, HA_COMM);
    sindx01=lj2/MAXSSIZE;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      if(SORD==1)MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_DOUBLE,group_size-1, 11*sindx01+sindx02, HA_COMM);
      else MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_FLOAT,group_size-1, 11*sindx01+sindx02, HA_COMM);
      sindx03+=MAXSSIZE;
    }
    if(SORD==1)MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,group_size-1, 11, HA_COMM);
    else MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,group_size-1, 11, HA_COMM);
    free(vecbivi);
    vecbivi=NULL;
    //vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
    //printf("rank %d color %d group11 %d\n",rank,color,group_size11);
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      MPI_Send(biviindx1+sindx03,MAXSSIZE, MPI_LONG,group_size-1, 12*sindx01+sindx02, HA_COMM);
      sindx03+=MAXSSIZE;
    }
    MPI_Send(biviindx1+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,group_size-1, 12, HA_COMM);
    free(obiviindx1);
    obiviindx1=NULL;
    //obiviindx1=realloc(obiviindx1,sizeof(long int));
    //printf("rank %d color %d group11 %d\n",rank,color,group_size11);
  }
  if(ha_id==group_size-1) {
    for(j6=0; j6<group_size-1; j6++) {
      //printf("!!!!!!!!!!!! rank11 %d nz0 %d j2 %d nz1 %d\n",rank,nz0,j2,nz1);
      MPI_Recv(&lj2,1, MPI_LONG,  MPI_ANY_SOURCE,10, HA_COMM,&status);
      i=status.MPI_SOURCE;
      obiviindx0=realloc(obiviindx0,(lj2+1)*sizeof(long int));
      obiviindx0[0]=-1;
      //int *obiviindx0= (int *) calloc (j2+1,sizeof(int));
      //obiviindx0[0]=-1;
      biviindx0=obiviindx0;
      biviindx0++;
      ha_cgetype *vecbivi0= (ha_cgetype *) calloc (lj2,sizeof(ha_cgetype));
      sindx01=lj2/MAXSSIZE;
      sindx03=0;
      //printf("!!!!!!!!!!!!1 rank11 %d nz0 %d j2 %d nz1 %d\n",rank,nz0,j2,nz1);
      for(sindx02=0; sindx02<sindx01; sindx02++) {
        if(SORD==1)MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_DOUBLE,i,11*sindx01+sindx02, HA_COMM,&status);
        else MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_FLOAT,i,11*sindx01+sindx02, HA_COMM,&status);
        //printf("rank %d sindx %d\n",rank,sindx02);
        sindx03+=MAXSSIZE;
      }
      if(SORD==1)MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,i,11, HA_COMM,&status);
      else MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,i,11, HA_COMM,&status);
      //printf("!!!!!!!!!!!!2 rank11 %d nz0 %d j2 %d nz1 %d\n",rank,nz0,j2,nz1);
      sindx03=0;
      for(sindx02=0; sindx02<sindx01; sindx02++) {
        MPI_Recv(biviindx0+sindx03,MAXSSIZE, MPI_LONG,i,12*sindx01+sindx02, HA_COMM,&status);
        sindx03+=MAXSSIZE;
      }
      MPI_Recv(biviindx0+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,i,12, HA_COMM,&status);
      //printf("!!!!!!!!!!!!2 rank11 %d nz0 %d j2 %d nz1 %d\n",rank,nz0,j2,nz1);
      my_spar_compl_(biviindx1,&nz0,biviindx0,&lj2,&nz1);
      vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
      obiviindx1=realloc(obiviindx1,(nz1+1)*sizeof(long int));
      obiviindx1[0]=-1;
      biviindx1=obiviindx1;
      biviindx1++;
      //printf("!!!!!!!!!!!!2 rank11 %d nz0 %d j2 %d nz1 %d\n",rank,nz0,j2,nz1);
      my_spar_add3l_(vecbivi,obiviindx1,&nz0,vecbivi0,obiviindx0,&lj2,&nz1);
      free(obiviindx0);
      obiviindx0=NULL;
      //obiviindx0=realloc(obiviindx0,sizeof(long int));
      //free(obiviindx0);
      free(vecbivi0);
      vecbivi0=NULL;
      printf("!!!!!!!!!!!! rank11 %d nz0 %ld j2 %ld nz1 %ld\n",rank,nz0,lj2,nz1);
      nz0=nz1;
    }
  }

  if(group_size11>1&&color==1) {
    //printf("OKKKK nz3 %d nz0 %d rank %d\n",nz3,nz0,rank);
    nz3=nz0;
    if(ha_id1!=group_size1-1) {
      lj2=nz3;
      //printf("OKKKK1 rank %d j2 %d\n",rank,j2);
      MPI_Send(&lj2,1, MPI_LONG,group_size1-1, 10, HA1_COMM);
      sindx01=lj2/MAXSSIZE;
      sindx03=0;
      for(sindx02=0; sindx02<sindx01; sindx02++) {
        if(SORD==1) MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_DOUBLE,group_size1-1, 11*sindx01+sindx02, HA1_COMM);
        else MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_FLOAT,group_size1-1, 11*sindx01+sindx02, HA1_COMM);
        sindx03+=MAXSSIZE;
      }
      if(SORD==1) MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,group_size1-1, 11, HA1_COMM);
      else MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,group_size1-1, 11, HA1_COMM);
      //printf("OKKKK1 rank %d j2 %d sz %ld\n",rank,j2,sizeof(obiviindx1));
      sindx03=0;
      for(sindx02=0; sindx02<sindx01; sindx02++) {
        MPI_Send(biviindx1+sindx03,MAXSSIZE, MPI_LONG,group_size1-1, 12*sindx01+sindx02, HA1_COMM);
        //printf("OKKKK1 rank %d s02 %d\n",rank,sindx02);
        sindx03+=MAXSSIZE;
      }
      //printf("OKKKK1 arnk %d nz1 %d %ld %ld %ld\n",rank,j2,obiviindx1[239829689],obiviindx1[239829689-1],obiviindx1[239829689+1]);
      MPI_Send(biviindx1+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,group_size1-1, 12, HA1_COMM);
      //printf("OKKKK1 rank %d\n",rank);
      free(obiviindx1);
      obiviindx1=NULL;
      //obiviindx1=realloc(obiviindx1,sizeof(long int));
      free(vecbivi);
      vecbivi=NULL;
      //vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
    }
    //printf("OKKKK rank %d\n",rank);
    if(ha_id1==group_size1-1) {
      for(j6=0; j6<group_size1-1; j6++) {
        //printf("OKKKK2 rank %d\n",rank);
        MPI_Recv(&lj2,1, MPI_LONG,  MPI_ANY_SOURCE,10, HA1_COMM,&status);
        i=status.MPI_SOURCE;
        //printf("OKKKK2 i %d rank %d j2 %d\n",i,rank,j2);
        obiviindx0=realloc(obiviindx0,(lj2+1)*sizeof(long int));
        obiviindx0[0]=-1;
        //int *obiviindx0= (int *) calloc (j2+1,sizeof(int));
        biviindx0=obiviindx0;
        biviindx0++;
        ha_cgetype *vecbivi0= (ha_cgetype *) calloc (lj2,sizeof(ha_cgetype));
        sindx01=lj2/MAXSSIZE;
        sindx03=0;
        //printf("OKKKK2 i %d rank %d s %d\n",i,rank,sindx01);
        for(sindx02=0; sindx02<sindx01; sindx02++) {
          if(SORD==1) MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_DOUBLE,i,11*sindx01+sindx02, HA1_COMM,&status);
          else MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_FLOAT,i,11*sindx01+sindx02, HA1_COMM,&status);
          //printf("rank %d sindx %d\n",rank,sindx02);
          sindx03+=MAXSSIZE;
        }
        if(SORD==1) MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,i,11, HA1_COMM,&status);
        else MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,i,11, HA1_COMM,&status);
        //printf("OKKKK2 i %d rank %d\n",i,rank);
        sindx03=0;
        for(sindx02=0; sindx02<sindx01; sindx02++) {
          MPI_Recv(biviindx0+sindx03,MAXSSIZE, MPI_LONG,i,12*sindx01+sindx02, HA1_COMM,&status);
          sindx03+=MAXSSIZE;
        }
        MPI_Recv(biviindx0+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,i,12, HA1_COMM,&status);
        //printf("OKKKK2 i %d rank %d nz0 %d nz1 %d\n",i,rank,nz0,j2);
        my_spar_compl_(biviindx1,&nz0,biviindx0,&lj2,&nz1);
        //printf("OKKKK2 i %d rank %d\n",i,rank);
        vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
        //printf("OKKKK2a i %d rank %d\n",i,rank);
        obiviindx1=realloc(obiviindx1,(nz1+1)*sizeof(long int));
        obiviindx1[0]=-1;
        //printf("OKKKK2b i %d rank %d\n",i,rank);
        biviindx1=obiviindx1;
        biviindx1++;
        //printf("OKKKK2 nz0 %d rank %d nz1 %d nz2 %d %ld %ld %ld\n",nz0,rank,j2,nz1,obiviindx0[239829689],obiviindx0[239829689-1],obiviindx0[239829689+1]);
        my_spar_add3l_(vecbivi,obiviindx1,&nz0,vecbivi0,obiviindx0,&lj2,&nz1);
        free(obiviindx0);
        obiviindx0=NULL;
        //obiviindx0=realloc(obiviindx0,sizeof(long int));
        //free(biviindx0);
        free(vecbivi0);
        printf("!!!!!!!!!!!! rank22 %d nz0 %ld j2 %ld nz1 %ld\n",rank,nz0,lj2,nz1);
        nz0=nz1;
      }
    }
  }
  /*if(rank!=mpisize-1)  {
    j2=nz3;
    //MPI_Bcast(&j2,1, MPI_INT,j6, PETSC_COMM_WORLD);
    MPI_Send(&j2,1, MPI_INT,mpisize-1, 10, PETSC_COMM_WORLD);
    //j3=j6;
    //MPI_Bcast(&j3,1, MPI_INT,j6, PETSC_COMM_WORLD);
    //if(rank==j6) {
      MPI_Send(vecbivi,j2, MPI_DOUBLE,mpisize-1, 11, PETSC_COMM_WORLD);
      MPI_Send(biviindx1,j2, MPI_INT,mpisize-1, 12, PETSC_COMM_WORLD);
      biviindx1=realloc(biviindx1,1*sizeof(int));
      vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
    //}
  }
    if(rank==mpisize-1) {
      for(j6=0; j6<mpisize-1; j6++){
      MPI_Recv(&j2,1, MPI_INT,  MPI_ANY_SOURCE,10, PETSC_COMM_WORLD,&status);
      i=status.MPI_SOURCE;
      //timestr1=clock();
      //printf("1 received %d numbers from 0. Message source = %d, " "tag = %d\n",1, i, status.MPI_TAG);
      //printf("j6 %d j2 %d\n",j6,j2);
      int *biviindx0= (int *) calloc (j2,sizeof(int));
      ha_cgetype *vecbivi0= (ha_cgetype *) calloc (j2,sizeof(ha_cgetype));
      MPI_Recv(vecbivi0,j2, MPI_DOUBLE,i,11, PETSC_COMM_WORLD,&status);
      //printf("j6 %d j2 %d\n",j6,j2);
      MPI_Recv(biviindx0,j2, MPI_INT,i,12, PETSC_COMM_WORLD,&status);
      //printf("j6 %d j2 %d nz0 %d\n",j6,j2,nz0);
      //nz2=nz0;
      my_spar_comp_(biviindx1,&nz0,biviindx0,&j2,&nz1);
      //printf("ok!!!!!!!!!!! nz1 %d\n",nz1);
      vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
      biviindx1=realloc(biviindx1,nz1*sizeof(int));
      my_spar_add_(vecbivi,biviindx1,&nz0,vecbivi0,biviindx0,&j2,&nz1);
      free(biviindx0);
      free(vecbivi0);
      printf("!!!!!!!!!!!! rank %d j %d j1 %d\n",rank,j,j1);
      nz0=nz1;
      //printf("Completed reduce Operation! Rank %d time %f\n",rank,((double)clock()-timestr1)/CLOCKS_PER_SEC);
      }
    }*/

  //printf("ok1!!!!!!!!!!! rank %d nz1 %d max %d\n",rank,nz1,name_len_max);
  //MPI_Barrier(PETSC_COMM_WORLD);
  /*for(j6=0; j6<mpisize-1; j6++) {
    j2=nz3;
    MPI_Bcast(&j2,1, MPI_INT,j6, PETSC_COMM_WORLD);
    j3=j6;
    MPI_Bcast(&j3,1, MPI_INT,j6, PETSC_COMM_WORLD);
    if(rank==j6) {
      MPI_Send(vecbivi,j2, MPI_DOUBLE,mpisize-1, 10, PETSC_COMM_WORLD);
      MPI_Send(biviindx1,j2, MPI_INT,mpisize-1, 11, PETSC_COMM_WORLD);
      biviindx1=realloc(biviindx1,1*sizeof(int));
      vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
    }
    //printf("j2 %d\n",j2);
    //printf("j1 %d rank000 %d\n",j6,rank);
    if(rank==mpisize-1) {
      int *biviindx0= (int *) calloc (j2,sizeof(int));
      ha_cgetype *vecbivi0= (ha_cgetype *) calloc (j2,sizeof(ha_cgetype));
      MPI_Recv(vecbivi0,j2, MPI_DOUBLE, j3, 10, PETSC_COMM_WORLD,&status);
      MPI_Recv(biviindx0,j2, MPI_INT, j3, 11, PETSC_COMM_WORLD,&status);
      nz2=nz0;
      my_spar_comp_(biviindx1,&nz2,biviindx0,&j2,&nz1);
      //if(biviindx1[nz2-1]<biviindx0[j2-1]) {
      //  j=0;
      //  //j1=0;
      //  i1=0;
      //  for(i=0; i<nz2; ++i) {
      //    while (biviindx1[i]>biviindx0[j])++j;//i1++;}
      //    //if(biviindx1[i]==biviindx0[j])j1++;
      //    if(biviindx1[i]<biviindx0[j])++i1;
      //  }
      //  nz1=j2+i1;//nz2+j2-j1;
      //  //printf("nz1 %d nz2 %d j1 %d j2 %d j %d i1 %d\n",nz1,nz2,j1,j2,j,i1);
      //} else {
      //  j=0;
      //  //j1=0;
      //  i1=0;
      //  for(i=0; i<j2; ++i) {
      //    while (biviindx0[i]>biviindx1[j]) ++j;//i1++;}
      //    //if(i<10)printf("i %d j %d x0 %d x1 %d",i,j,biviindx0[i],biviindx1[j]);
      //    if(biviindx0[i]<biviindx1[j])++i1;
      //    //if(biviindx0[i]==biviindx1[j])j1++;
      //  }
      //  nz1=nz2+i1;//nz2+j2-j1;
      //}
      vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
      biviindx1=realloc(biviindx1,nz1*sizeof(int));
      //timestr=clock();
      my_spar_add_(vecbivi,biviindx1,&nz2,vecbivi0,biviindx0,&j2,&nz1);
  //      j=j2-1;
  //      j1=nz2-1;
  //      i=nz1-1;
  //      while(i>-1) {
  //        if(biviindx1[j1]==biviindx0[j]) {
  //          vecbivi[i]=vecbivi[j1]+vecbivi0[j];
  //          biviindx1[i]=biviindx1[j1];
  //          --j1;
  //          --j;
  //          --i;
  //          continue;
  //        }
  //        i1=i;
  //        while(biviindx1[j1]>biviindx0[j]) {
  //          //vecbivi[i]=vecbivi[j1];
  //          //biviindx1[i]=biviindx1[j1];
  //          //i1;
  //          --j1;
  //          --i;
  //          if(j1<=-1)break;
  //        }
  //        if(i<i1){
  //          memmove(biviindx1+i+1,biviindx1+j1+1,(i1-i)*sizeof(unsigned int));
  //          memmove(vecbivi+i+1,vecbivi+j1+1,(i1-i)*sizeof(ha_cgetype));
  //          //printf("i1 %d i %d indx1 %d indx %d\n",i1,i,biviindx1[i+1+i-i1],biviindx1[i1+1+i-i1]);
  //        }
  //        i1=i;
  //        while(biviindx1[j1]<biviindx0[j]) {
  //          //vecbivi[i]=vecbivi0[j];
  //          //biviindx1[i]=biviindx0[j];
  //          --j;
  //          --i;
  //          if(j<=-1)break;
  //        }
  //        if(i<i1){
  //          memcpy(biviindx1+i+1,biviindx0+j+1,(i1-i)*sizeof(unsigned int));
  //          memcpy(vecbivi+i+1,vecbivi0+j+1,(i1-i)*sizeof(ha_cgetype));
  //          //printf("i1 %d i %d indx1 %d indx %d\n",i1,i,biviindx1[i+1+i-i1],biviindx1[i1+1+i-i1]);
  //        }
  //      }
      //printf("Vec add! Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
      free(biviindx0);
      free(vecbivi0);
      printf("!!!!!!!!!!!!i %d j %d j1 %d\n",i,j,j1);
      nz0=nz1;
    }
  }*/
  //free(biviindx0);
  //free(vecbivi0);
  //if(rank!=mpisize-1) {
  //vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
  //biviindx1=realloc(biviindx1,1*sizeof(unsigned int));
  //}
  forint vecbiuisize=0;
  vecbiuisize=VecSize-sumrowcolin;
  printf("after bivi rank %d\n",rank);
  HaReduceNoComp(vecbiui,vecbiuisize,mpisize,rank,mpisize-1);
  printf("Completed MPI_Reduce Operation! Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  if(rank!=mpisize-1){
    free(vecbiui);
    vecbiui=NULL;
    //vecbiui=realloc(vecbiui,sizeof(ha_cgetype));//free(vecbiui);
  }
  timestr=clock();
  ha_cgetype *xd;//= (ha_cgetype*)calloc(vecbiuisize,sizeof(ha_cgetype));
  ha_cgetype *x0;//,*vecbivi0;//= (ha_cgetype*)calloc(VecSize,sizeof(ha_cgetype));
  long int j7,lnz;
  if(rank==mpisize-1) {
    //for(i=0; i<vecbiuisize; i++)yd[i]-=vecbiui[i];
    //vecbiui=realloc(vecbiui,sizeof(ha_cgetype));//free(vecbiui);
    Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatD[0]->data;//*aa=subA->data;
    ai= aa->i;
    aj= aa->j;
    vals=aa->a;
    lnz=aa->nz;
    nrow=submatD[0]->rmap->n;
    ncol=submatD[0]->cmap->n;
    printf("D ncol %d nrow %d nz %ld\n",ncol,nrow,lnz);
    obiviindx0=realloc(obiviindx0,(lnz+1)*sizeof(long int));
    obiviindx0[0]=-1;
    //obiviindx0[0]=-1;
    biviindx0=obiviindx0;
    biviindx0++;
    //int *biviindx0=(int *) calloc (nz,sizeof(int));
    lj2=0;
    for(li=0; li<nrow-1; li++) {
      j7=li*ncol;
      for(j=ai[li]; j<ai[li+1]; j++) {
        biviindx0[lj2]=aj[j]+j7;
        //if(biviindx0[j2]<0)printf("biviindx0[j2] %ld j %d j3 %d\n",biviindx0[j2],j,j3);
        lj2++;
      }
    }
    j7=li*ncol;
    for(j=ai[li]; j<lnz; j++) {
      biviindx0[lj2]=aj[j]+j7;
      //if(biviindx0[j2]<0)printf("biviindx0[j2] %ld j %d j3 %d\n",biviindx0[j2],j,j3);
      lj2++;
    }
    //j=nz-1;
    //j1=0;
    my_spar_compl_(biviindx1,&nz0,biviindx0,&lnz,&nz1);
    printf("nz0 %ld nz %ld nz1 %ld\n",nz0,lnz,nz1);
    /*if(biviindx1[nz0-1]<biviindx0[j2-1]) {
      j=0;
      j1=0;
      for(i=0; i<nz0; i++) {
        while (biviindx1[i]>biviindx0[j])j++;
        //if(biviindx1[i]==biviindx0[j])j1++;
        if(biviindx1[i]<biviindx0[j])j1++;
      }
      nz1=j2+j1;//nz0+j2-j1;
      //printf("nz1 %d j1 %d j2 %d j %d\n",nz1,j1,j2,j);
    } else {
      j=0;
      j1=0;
      for(i=0; i<j2; i++) {
        while (biviindx0[i]>biviindx1[j]) j++;
        //if(i<10)printf("i %d j %d x0 %d x1 %d",i,j,biviindx0[i],biviindx1[j]);
        //if(biviindx0[i]==biviindx1[j])j1++;
        if(biviindx0[i]<biviindx1[j])j1++;
      }
      nz1=nz0+j1;//nz0+j2-j1;
      //printf("nz1 %d j1 %d nz2 %d j %d j2 %d i %d\n",nz1,j1,nz2,j,j2,i);
    }*/



    /*for(i=nz0-1; i>-1; i--) {
      //printf("i %d j1 %d j %d v %d v0 %d v1 %d\n",i,j1,j,biviindx[i],biviindx0[j],biviindx[j1+1]);
      if(biviindx1[i]>biviindx0[j]) {
        j1--;
        continue;
      }
      if(biviindx1[i]==biviindx0[j]) {
        j--;
        j1--;
        if(j==-1) {
          i--;
          break;
        }
        continue;
      }
      for(j3=j; j3>-1; j3--) {
        if(biviindx1[i]>biviindx0[j3]) {
          j1--;
          break;
        }
        if(biviindx1[i]==biviindx0[j3]) {
          j3--;
          j1--;
          break;
        }
        j1--;
      }
      j=j3;
      if(j==-1)break;
    }
    //printf("j1 %d j %d i %d v %d v0 %d\n",j1,j,i,biviindx1[i+1],biviindx0[j+1]);
    if(j>-1)j1-=j+1;//for(i=j;i>-1;i--){
    if(i>-1)j1-=i+1;//break
    nz1=-j1;*/
    printf("nz1 %ld\n",nz1);
//    if(!iter) {
    ldsize=ceil((laD/100.0)*nz1);
    vecbivi=realloc(vecbivi,ldsize*sizeof(ha_cgetype));
    //biviindx1=realloc(biviindx1,nz1*sizeof(unsigned int));
    int *irn1=(int *) calloc (nz1,sizeof(int));
    int *jcn=(int *) calloc (ldsize,sizeof(int));
    printf("indx0 %ld indx1 %ld\n",obiviindx1[0],obiviindx0[0]);
    my_spar_add4l_(vecbivi,obiviindx1,irn1,jcn,&nz0,vals,obiviindx0,&lnz,&nz1,&ncol);//my_spar_add1l_
    /*j=nz-1;
    j1=nz0-1;
    i=nz1-1;
    while(i>-1) {
      if(biviindx1[j1]==biviindx0[j]) {
        vecbivi[i]=vecbivi[j1]+vals[j];
        irn1[i]=biviindx1[j1]/ncol+1;
        jcn[i]=biviindx1[j1]%ncol+1;//biviindx1[j1]-(irn1[i]-1)*ncol+1;
        --j1;
        --j;
        --i;
        continue;
      }
      i1=i;
      while(biviindx1[j1]>biviindx0[j]) {
        //vecbivi[i]=vecbivi[j1];
        irn1[i]=biviindx1[j1]/ncol+1;
        jcn[i]=biviindx1[j1]%ncol+1;//biviindx1[j1]-irn1[i]*ncol+1;
        //irn1[i]++;
        //i1;
        --j1;
        --i;
        if(j1<=-1)break;
      }
      if(i<i1){
        memmove(vecbivi+i+1,vecbivi+j1+1,(i1-i)*sizeof(ha_cgetype));
      }
      i1=i;
      while(biviindx1[j1]<biviindx0[j]) {
        //vecbivi[i]=vals[j];
        irn1[i]=biviindx0[j]/ncol+1;
        jcn[i]=biviindx0[j]%ncol+1;//biviindx0[j]-(irn1[i]-1)*ncol+1;
        --j;
        --i;
        if(j<=-1)break;
      }
      if(i<i1){
        memcpy(vecbivi+i+1,vals+j+1,(i1-i)*sizeof(ha_cgetype));
      }
    }*/

    /*j=nz-1;
    j1=nz0-1;
    for(i=nz1-1; i>-1; i--) {
      if(biviindx1[j1]>biviindx0[j]) {
        vecbivi[i]=vecbivi[j1];
        //biviindx1[i]=biviindx1[j1];
        irn1[i]=biviindx1[j1]/ncol;
        jcn[i]=biviindx1[j1]-irn1[i]*ncol+1;
        irn1[i]++;
        j1--;
      } else {
        if(biviindx1[j1]<biviindx0[j]) {
          vecbivi[i]=vals[j];
          //biviindx1[i]=biviindx0[j];
          irn1[i]=biviindx0[j]/ncol+1;
          jcn[i]=biviindx0[j]-(irn1[i]-1)*ncol+1;
          j--;
        } else {
          vecbivi[i]=vecbivi[j1]+vals[j];
          //biviindx1[i]=biviindx1[j1];
          irn1[i]=biviindx1[j1]/ncol+1;
          jcn[i]=biviindx1[j1]-(irn1[i]-1)*ncol+1;
          j1--;
          j--;
        }
      }
    }*/


    /*j=nz-1;
    j1=nz1-1;
    printf("anz2 %d nz %d\n",nz0,nz);
    for(i=nz0-1; i>-1; i--) {
      if(biviindx1[i]>biviindx0[j]) {
        vecbivi[j1]=vecbivi[i];
        biviindx1[j1]=biviindx1[i];
        j1--;
        continue;
      }
      if(biviindx1[i]==biviindx0[j]) {
        vecbivi[j1]=vecbivi[i]+vals[j];
        biviindx1[j1]=biviindx1[i];
        j--;
        j1--;
        if(j==-1)break;
        continue;
      }
      for(j3=j; j3>-1; j3--) {
        if(biviindx1[i]>biviindx0[j3]) {
          vecbivi[j1]=vecbivi[i];
          biviindx1[j1]=biviindx1[i];
          j1--;
          break;
        }
        if(biviindx1[i]==biviindx0[j3]) {
          vecbivi[j1]=vecbivi[i]+vals[j3];
          biviindx1[j1]=biviindx1[i];
          j3--;
          j1--;
          break;
        }
        biviindx1[j1]=biviindx0[j3];
        vecbivi[j1]=vals[j3];
        j1--;
      }
      j=j3;
      if(j==-1)break;
    }
    if(j>-1)for(i=j; i>-1; i--) {
        biviindx1[j1]=biviindx0[i];
        vecbivi[j1]=vals[i];
        j1--;
      }*/
    //nz0=nz1;
    MatDestroy(&submatD[0]);
    //free(biviindx0);
    lnz=nz1;
    //int *irn1=(int *) calloc (laD*nz,sizeof(int));
    //int *jcn=(int *) calloc (laD*nz,sizeof(int));
    //int *nnz=(int *) calloc (nrow,sizeof(int));
    //for(i=0; i<nrow; i++)nnz[i]=1;
    /*for(i=0; i<nz; i++) {
      j3=biviindx1[i]/ncol;
      //nnz[j3]++;
      //printf("nrow %d ncol %d j3 %d nnz %d",nrow,ncol,j3,nnz[j3]);
      irn1[i]=j3+1;
      jcn[i]=biviindx1[i]-j3*ncol+1;
    }*/
    //obiviindx1=realloc(obiviindx1,sizeof(long int));
    free(obiviindx1);
    obiviindx1=NULL;
    //obiviindx0=realloc(obiviindx0,sizeof(long int));
    free(obiviindx0);
    obiviindx0=NULL;
    irn1=realloc(irn1,ldsize*sizeof(int));
    int *insizeD=(int *) calloc (6,sizeof(int));
    insizeD[0]=nrow;
    insizeD[1]=ncol;
    insizeD[2]=lnz;
    insizeD[3]=laD;
    insizeD[5]=ldsize;
    printf("ok\n");
    xd=(ha_cgetype *) calloc (vecbiuisize,sizeof(ha_cgetype));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
    printf("ok!!!!!!!! vecbiuisize %ld nrow %d\n",vecbiuisize,nrow);
    spec48_ssol2la_(insizeD,irn1,jcn,vecbivi,vecbiui,xd);
    printf("ok!!!!!!!!\n");


    /*//Beg petsc
    int *nnz=(int *) calloc (nrow,sizeof(int));
    for(i=0; i<nrow; i++)nnz[i]=1;
    for(i=0; i<nz; i++) nnz[irn1[i]-1]++;
    Vec bright,b1,xsol;
    Mat      A,C,D,F;
    KSP      ksp;
    PC       pc;
    PetscInt  nzperrow,its;
    PetscErrorCode ierr;
    PetscScalar *array;
    nzperrow=0;
    for(i=0; i<nrow; i++)if(nzperrow<nnz[i])nzperrow=nnz[i];
    //for(i=0; i<nrow; i++)
    printf("nzperrow %d\n",nzperrow);
    ierr = VecCreate(PETSC_COMM_SELF,&bright);
    VecSetType(bright,VECSEQ);
    VecSetSizes(bright,ncol,ncol);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(bright);
    CHKERRQ(ierr);
    //ierr = PetscMalloc(nrow*sizeof(PetscInt),&nnz);
    //CHKERRQ(ierr);
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nrow,ncol,nzperrow,nnz,&A);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);
    CHKERRQ(ierr);
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nrow,ncol,nzperrow,nnz,&C);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(C);
    CHKERRQ(ierr);
    //ierr = PetscFree(nnz);
    //CHKERRQ(ierr);
    free(nnz);
    for(i=0; i<nz; i++){
      //printf("irn1 %d\n",irn1[i]);
      //irn1[i]--;
      //printf("irn1 %d\n",irn1[i]);
      //jcn[i]--;
      if(vecbivi[i]>1e-2||vecbivi[i]<-1e-2)ierr = MatSetValue(C,irn1[i]-1,jcn[i]-1,(PetscScalar)vecbivi[i],INSERT_VALUES);
      ierr = MatSetValue(A,irn1[i]-1,jcn[i]-1,(PetscScalar)vecbivi[i],INSERT_VALUES);
      //printf("irn1 %d jcn %d v %f\n",irn1[i],jcn[i],vecbivi[i]);
    }
    //for(i=0; i<10; i++)printf("irn %d jcn %d A %lf\n",irn1[i],jcn[i],vecbivi[i]);
    for(i=0; i<nrow; i++)
      ierr = MatSetValue(A,i,i,0.0,ADD_VALUES);
    for(i=0; i<nrow; i++)
      ierr = MatSetValue(C,i,i,0.0,ADD_VALUES);
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    printf("nzperrow %d\n",nzperrow);
    for(i=0; i<nrow; i++)VecSetValues(bright,1,&i,&vecbiui[i],INSERT_VALUES);//printf("yd %f\n",yd[i]);}
    printf("nzperrow %d\n",nzperrow);
    ierr = VecAssemblyBegin(bright);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(bright);
    CHKERRQ(ierr);
    ierr = VecDuplicate(bright,&xsol);
    CHKERRQ(ierr);
    //ierr = VecDuplicate(bright,&b1);
    //CHKERRQ(ierr);
    //MatView(A,0);
    //VecView(bright,0);
    PetscScalar one=1.0;
    ierr = VecSet(xsol,one);
    CHKERRQ(ierr);
    //spec48_ssol2la_(insizeD,irn1,jcn,vecbivi,vecbiui,xd);
    //for(i=0; i<nrow; i++)VecSetValues(xsol,1,&i,&xd[i],INSERT_VALUES);
    //ierr = VecAssemblyBegin(xsol);
    //CHKERRQ(ierr);
    //ierr = VecAssemblyEnd(xsol);
    //CHKERRQ(ierr);
    //printf("OKKK!!!\n");
    //VecView(xsol,0);
    //MatMult(A,xsol,bright);
    //VecView(bright,0);
    //ierr = MatMultTranspose(A,bright,b1);
    printf("OKKK!!!\n");
    //ierr = MatTranspose(A,MAT_INITIAL_MATRIX,&C);
    //ierr = MatMatMult(C,A,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&D);
    //VecView(b1,0);

    ierr = KSPCreate(PETSC_COMM_SELF,&ksp);
    CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,A,C,SAME_PRECONDITIONER);
    CHKERRQ(ierr);
    //KSPSetType(ksp,KSPPREONLY);
    KSPSetType(ksp,KSPBCGS);
    //KSPSetType(ksp,KSPFGMRES);
    ierr=KSPGetPC(ksp,&pc);
    CHKERRQ(ierr);
    ierr=PCSetType(pc,PCLU);
    //ierr=PCSetType(pc,PCJACOBI);
    CHKERRQ(ierr);
    //PCFactorSetMatOrderingType(pc,MATORDERINGND);
    //PCFactorSetShiftType(pc, MAT_SHIFT_POSITIVE_DEFINITE);
    //PCFactorSetAllowDiagonalFill(pc);
    //PCFactorReorderForNonzeroDiagonal(pc,1.e-5);
    //PCFactorSetShiftAmount(pc,PETSC_DECIDE);
    //KSPSetType(ksp,KSPCGNE);
    //PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU);
    //PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU_DIST);
    PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
    PCFactorSetUpMatSolverPackage(pc);
    ierr = PCFactorGetMatrix(pc,&F);CHKERRQ(ierr);
    //ierr = MatMumpsSetIcntl(F,6,2);CHKERRQ(ierr);
    ierr = MatMumpsSetIcntl(F,7,5);CHKERRQ(ierr);
    CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,1e-10/nrow,1e-50,PETSC_DEFAULT,PETSC_DEFAULT);
    CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);
    printf("OKKK!!!\n");
    ierr = KSPSolve(ksp,bright,xsol);
    CHKERRQ(ierr);
    printf("ok!!!!!!!!\n");
    ierr = KSPGetIterationNumber(ksp,&its);
    CHKERRQ(ierr);
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(ksp,&reason);
    CHKERRQ(ierr);
    printf("Reason %d \n",reason);
    if (reason==KSP_DIVERGED_INDEFINITE_PC) {
      printf("\nDivergence because of indefinite preconditioner;\n");
      printf("Run the executable again but with '-pc_factor_shift_type POSITIVE_DEFINITE' option.\n");
    } else if (reason<0) {
      printf("\nOther kind of divergence: this should not happen.\n");
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\nConvergence in %d iterations.\n",its);
    }
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF,"Iterations %d\n",its);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);
    //ierr = MatDestroy(&F);
    //CHKERRQ(ierr);
    ierr = VecDestroy(&bright);
    CHKERRQ(ierr);
    printf("Reason %d \n",reason);
    //PCDestroy(&pc);
    ierr = KSPDestroy(&ksp);
    CHKERRQ(ierr);
    printf("Reason %d \n",reason);
    VecGetArray(xsol,&array);
    memcpy (xd,array,nrow*sizeof(PetscScalar) );
    VecRestoreArray(xsol,&array);
    ierr = VecDestroy(&xsol);
    CHKERRQ(ierr);
    //End petsc*/

    free(vecbivi);
    vecbivi=NULL;
    //vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
    free(vecbiui);
    vecbiui=NULL;
    //vecbiui=realloc(vecbiui,1*sizeof(ha_cgetype));
    //free(yd);
    free(irn1);
    free(jcn);
    //free(values);
    free(insizeD);
//    }
//    else {
//      vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
//      //biviindx1=realloc(biviindx1,nz1*sizeof(unsigned int));
//      int *irn1=(int *) calloc (nz1,sizeof(int));
//      int *jcn=(int *) calloc (nz1,sizeof(int));
//      ha_cgetype *vecbivi1=(ha_cgetype *) calloc (nz1,sizeof(ha_cgetype));
//      int *irn2=(int *) calloc (nz1,sizeof(int));
//      int *jcn1=(int *) calloc (nz1,sizeof(int));
//      my_spar_add2_(vecbivi,obiviindx1,irn1,jcn,&nz0,vals,obiviindx0,&nz,&nz1,&ncol,vecbivi1,irn2,jcn1,&j2int,&cntl3);
//      j2=j2int;
//
//      /*j=nz-1;
//      j1=nz0-1;
//      i=nz1-1;
//      j2=nz1-1;
//      while(i>-1) {
//        if(biviindx1[j1]==biviindx0[j]) {
//          vecbivi[i]=vecbivi[j1]+vals[j];
//          irn1[i]=biviindx1[j1]/ncol+1;
//          jcn[i]=biviindx1[j1]%ncol+1;//biviindx1[j1]-(irn1[i]-1)*ncol+1;
//          if(vecbivi[i]>cntl3||vecbivi[i]<-cntl3){
//            vecbivi1[j2]=vecbivi[i];
//            irn2[j2]=irn1[i];
//            jcn1[j2]=jcn[i];
//            j2--;
//          }
//          --j1;
//          --j;
//          --i;
//          continue;
//        }
//        i1=i;
//        while(biviindx1[j1]>biviindx0[j]) {
//          //vecbivi[i]=vecbivi[j1];
//          irn1[i]=biviindx1[j1]/ncol+1;
//          jcn[i]=biviindx1[j1]%ncol+1;//biviindx1[j1]-irn1[i]*ncol+1;
//          //irn1[i]++;
//          if(vecbivi[j1]>cntl3||vecbivi[j1]<-cntl3){
//            vecbivi1[j2]=vecbivi[j1];
//            irn2[j2]=irn1[i];
//            jcn1[j2]=jcn[i];
//            j2--;
//          }
//          --j1;
//          --i;
//          if(j1<=-1)break;
//        }
//        if(i<i1){
//          memmove(vecbivi+i+1,vecbivi+j1+1,(i1-i)*sizeof(ha_cgetype));
//        }
//        i1=i;
//        while(biviindx1[j1]<biviindx0[j]) {
//          //vecbivi[i]=vals[j];
//          irn1[i]=biviindx0[j]/ncol+1;
//          jcn[i]=biviindx0[j]%ncol+1;//biviindx0[j]-(irn1[i]-1)*ncol+1;
//          if(vals[j]>cntl3||vals[j]<-cntl3){
//            vecbivi1[j2]=vals[j];
//            irn2[j2]=irn1[i];
//            jcn1[j2]=jcn[i];
//            j2--;
//          }
//          --j;
//          --i;
//          if(j<=-1)break;
//        }
//        if(i<i1){
//          memcpy(vecbivi+i+1,vals+j+1,(i1-i)*sizeof(ha_cgetype));
//        }
//      }*/
//      nz2=nz1-j2-1;
//      obiviindx1=realloc(obiviindx1,sizeof(int));
//      memmove(vecbivi1,vecbivi1+j2+1,nz2*sizeof(ha_cgetype));
//      memmove(irn2,irn2+j2+1,nz2*sizeof(int));
//      memmove(jcn1,jcn1+j2+1,nz2*sizeof(int));
//      vecbivi1=realloc(vecbivi1,nz2*sizeof(ha_cgetype));
//      irn2=realloc(irn2,nz2*sizeof(int));
//      jcn1=realloc(jcn1,nz2*sizeof(int));
//      vecbivi1=realloc(vecbivi1,nz2*laD*sizeof(ha_cgetype));
//      irn2=realloc(irn2,nz2*laD*sizeof(int));
//      jcn1=realloc(jcn1,nz2*laD*sizeof(int));
//      /*j=1;
//      for(i=1;i<nz1;i++){
//        if(irn1[i-1]!=irn1[i]){irn1[j]=i+1;j++;}
//      }
//      if(j!=nrow)printf("Matrix is singular!\n");
//      else {
//        irn1[nrow]=nz1+1;
//      }
//      irn1=realloc(irn1,(nrow+1)*sizeof(int));*/
//      MatDestroy(&submatD[0]);
//      obiviindx0=realloc(obiviindx0,sizeof(int));
//      //free(biviindx0);
//      nz=nz1;
//      int *insizeD=(int *) calloc (3,sizeof(int));
//      insizeD[0]=nrow;
//      insizeD[1]=nz;
//      insizeD[2]=nz2;
//      insizeD[3]=laD;
//      printf("ok vecbivi %lf\n",vecbivi[nz1-1]);
//      xd=(ha_cgetype *) calloc (vecbiuisize,sizeof(ha_cgetype));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
//      printf("ok!!!!!!!! vecbiuisize %ld nrow %d\n",vecbiuisize,nrow);
//      my_iter_mi24ds_(insizeD,irn1,jcn,vecbivi,irn2,jcn1,vecbivi1,vecbiui,xd,&cntl3);
//      //for(i=0;i<nrow;i++)printf("xd %lf\n",xd[i]);
//      printf("ok!!!!!!!!\n");
//      vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
//      vecbiui=realloc(vecbiui,1*sizeof(ha_cgetype));
//      free(irn1);
//      free(jcn);
//      free(vecbivi1);
//      free(irn2);
//      free(jcn1);
//      free(insizeD);
//    }
    x0=(ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));//realloc (x0,*sizeof(ha_cgetype));
    for(i=0; i<vecbiuisize; i++) {
      //VecGetValues(xsol,1,&i,&xd[i]);
      x0[indicesC[i]]=xd[i];
      //printf("i %d x1d %lf\n",i,xd[i]);
    }
    //VecView(xsol,0);
    //ierr = VecDestroy(&xsol);
    //CHKERRQ(ierr);
  }
  free(vecbiui);
    obiviindx1=realloc(obiviindx1,sizeof(long int));
    free(obiviindx1);
    obiviindx1=NULL;
    obiviindx0=realloc(obiviindx0,sizeof(long int));
    free(obiviindx0);
    obiviindx0=NULL;
  //free(obiviindx1);
  //free(obiviindx0);
  free(vecbivi);
  free(indicesC);
  if(rank!=mpisize-1) {
    //free(bivinzrow);
    //free(bivinzcol);
    //free(bivinzrow1);
    //free(bivinzcol1);
    xd=(ha_cgetype *) calloc (vecbiuisize,sizeof(ha_cgetype));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
    x0=(ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));//realloc (x0,*sizeof(ha_cgetype));
  }
  PetscFree(submatD);
  //PetscFree(submatCT);
  MPI_Barrier(PETSC_COMM_WORLD);
  if(SORD==1)MPI_Bcast(xd, vecbiuisize, MPI_DOUBLE,mpisize-1, PETSC_COMM_WORLD);
  else MPI_Bcast(xd, vecbiuisize, MPI_FLOAT,mpisize-1, PETSC_COMM_WORLD);
  printf("Interface problem solution Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  //printf("sizeof forint %d sizeof int %d sizeof PesctInt %d sof PetscScalar %d sof ha_cgetype %d\n",sizeof(forint),sizeof(int),sizeof(PetscInt),sizeof(PetscScalar),sizeof(ha_cgetype));
  timestr=clock();
  xi1indx=0;
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      ha_cgetype *biui0= (ha_cgetype *) calloc (ha_ndblocks[j1+begblock[rank]],sizeof(ha_cgetype));
//      if(j1<0) {
//        strcpy(tempchar,"unzip ");
//        strcat(tempchar,fn01[j1]);
//        strcat(tempchar,".zip");
//        system(tempchar);
//        strcpy(tempchar,fn01[j1]);
//        strcat(tempchar,".zip");
//        remove(tempchar);
//        strcpy(tempchar,"unzip ");
//        strcat(tempchar,fn02[j1]);
//        strcat(tempchar,".zip");
//        system(tempchar);
//        strcpy(tempchar,fn02[j1]);
//        strcat(tempchar,".zip");
//        remove(tempchar);
//        strcpy(tempchar,"unzip ");
//        strcat(tempchar,fn03[j1]);
//        strcat(tempchar,".zip");
//        system(tempchar);
//        strcpy(tempchar,fn03[j1]);
//        strcat(tempchar,".zip");
//        remove(tempchar);
//      }

      fp1 = fopen(fn01[j1], "rb");
      //printf("fn1 %s\n",fn01[j1]);
      if (fp1==NULL)printf("File opening error\n");
      fp2 = fopen(fn02[j1], "rb");
      if (fp2==NULL)printf("File opening error\n");
      fp3 = fopen(fn03[j1], "rb");
      if (fp3==NULL)printf("File opening error\n");
      la1=ceil((insize[j1*insizes+9]/100.0)*insize[j1*insizes+2]);
      int *irne = (int*)calloc(la1,sizeof(int));
      int *keep = (int*)calloc(insize[j1*insizes+12],sizeof(int));
      ha_cgetype *vale = (ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
      freadresult=fread(irne,sizeof(int),la1,fp1);
      freadresult=fread(keep,sizeof(int),insize[j1*insizes+12],fp2);
      freadresult=fread(vale,sizeof(ha_cgetype),la1,fp3);
      fclose(fp1);
      fclose(fp2);
      fclose(fp3);
      remove(fn01[j1]);
      remove(fn02[j1]);
      remove(fn03[j1]);
      //for(i=0;i<la;i++)printf("vav %f\n",vale[i]);
      Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatC[j1]->data;//*aa=subA->data;
      //printf("0 rank %d j1 %d vec %lf\n",rank,j1,vecbivi[22800]);
      ai= ac->i;
      aj= ac->j;
      vals=ac->a;
      nz=ac->nz;
      nrow=submatC[j1]->rmap->n;
      ncol=submatC[j1]->cmap->n;
      //printf("nrow %d\n",nrow);
      ha_cgetype *be0 = (ha_cgetype*)calloc(nrow,sizeof(ha_cgetype));
      for(i=0; i<nrow-1; i++) {
        for(j=ai[i]; j<ai[i+1]; j++) {
          be0[i]+=vals[j]*xd[aj[j]];
        }
      }
      for(j=ai[i]; j<nz; j++) {
        be0[i]+=vals[j]*xd[aj[j]];
      }
      MatDestroy(&submatC[j1]);
      if(insize[j1*insizes+16]!=la1)insize[j1*insizes+16]=la1;
      spec48m_esol_(insize+j1*insizes,irne,vale,keep,be0,biui0);
      //printf("here!\n");
      free(irne);
      free(keep);
      free(vale);
      free(be0);
      ISGetIndices(colindices[j1],&nindices);
      xi1point=xi1+xi1indx;
      for(i=0; i<ha_ndblocks[j1+begblock[rank]]; i++) {
        x0[nindices[i]]+=xi1point[i]-biui0[i];//vals[i]-biui0[i];//vecval;
      }
      ISRestoreIndices(colindices[j1],&nindices);
      free(biui0);
      xi1indx+=ha_ndblocks[j1+begblock[rank]];
    }
    //MPI_Barrier(PETSC_COMM_WORLD);
  }
  free(xi1);
  free(insize);
  ierr = PetscFree(submatC);
  //free(vnz);
  free(begblock);
  printf("Solution calculation Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  //printf("Ok0 rank %d\n",rank);
  MPI_Barrier(PETSC_COMM_WORLD);
  //printf("Ok0 rank %d\n",rank);
  free(xd);
  timestr=clock();
  if(SORD==1)MPI_Allreduce(x0,x1,VecSize, MPI_DOUBLE, MPI_SUM,PETSC_COMM_WORLD);
  else MPI_Allreduce(x0,x1,VecSize, MPI_FLOAT, MPI_SUM,PETSC_COMM_WORLD);
  printf("Reduce solution Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  free(x0);
  for (i=0; i<nmatinplus; i++){
    free(fn01[i]);
    free(fn02[i]);
    free(fn03[i]);
  }
  free(fn01);
  free(fn02);
  free(fn03);
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&colindices[i]);
    CHKERRQ(ierr);
  }
  PetscFree(colindices);
  printf("OK rui rank %d\n",rank);
  return 0;
}

int HaNDBBDParPre(Mat A, Vec b, ha_cgetype *x1, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, uvadd ndblock,uvadd nreg,uvadd ntime,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,uvdim laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol) {//,bool iter
  IS *rowindices=NULL,*colindices=NULL,*colindicesbc1=NULL,*rowBBij=NULL,*colBBij=NULL;//,*colindicesbcpm,*colindicesbcpm1
  //const PetscInt *nindices;
  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nmatint,nmatinplust,nmatminust,nrowcolin,i,i1,j,j0,j1,j2,j3,j4,j5,j6,j7,l0,l1,l2,l3,l4,l5,rank1,proc1=0,nnzmax,j1nz,j1irnbs;//,sumrowcolin
  Mat *submatAij=NULL,*submatBBij=NULL;//,*submatCij,*submatBij;,*submatB
  Vec vecxd;
  PetscInt *ai,*aj,*a1i,*a1j,*a2i,*a2j,*aic,*ajc;
  PetscReal cntl6in;
  uvadd lasize,ldsize;
  //long int *biviindx1,*biviindx0,j7l;
  //long int *obiviindx1= (long int*)calloc(1,sizeof(long int));
  //long int *obiviindx0= (long int*)calloc(1,sizeof(long int));
  size_t freadresult,fwrt;
  //obiviindx1[0]=-1;
  //obiviindx0[0]=-1;
  PetscInt nrow,ncol,nz,nrowc,ncolc,nrowb,ncolb,nzc,nzv;
  MPI_Status   status;
  //clock_t timestr,timeend;
  int j2int,la1,fd1,fd2,fd3,unequal=0;
  bool ifremove=false;
  PetscScalar *vals,*val1s,*val2s,*valsc,vecval;
  ha_cgetype vval;
  PetscErrorCode ierr;
  PetscViewer viewer;
  FILE* presolfile;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  int *ndblockinrank= (int *) calloc (mpisize,sizeof(int));
  char filename[1024],rankname[1024],j1name[1024],tempchar[1024];//,fn1[1024],fn2[1024],fn3[1024]
  if(rank<10)strcpy(rankname,"000");
  if(rank<100&&rank>9)strcpy(rankname,"00");
  if(rank<1000&&rank>99)strcpy(rankname,"0");
  if(rank>=1000)rankname[0]='\0';
  sprintf(filename, "%d",rank);
  strcat(rankname,filename);
  int nmatinBB=1;
  nmatint=(uvadd)ntime/mpisize;
  if(nmatint*mpisize<ntime)unequal=1;
  nmatminust=nmatint;
  if(nmatint*mpisize<ntime)nmatinplust=nmatint+1;
  else nmatinplust=nmatint;
  for(i=0; i<mpisize; i++)if(rank+1<=ntime-mpisize*nmatint)nmatint++;
  nmatin=(nreg+1)*nmatint;
  nmatinplus=(nreg+1)*nmatinplust;
  begblock[rank]=nmatin;
  printf("rank %d nmatin %d nmatint %d nmplus %d\n",rank,nmatin,nmatint,nmatinplus);
  for(i=0; i<mpisize; i++) {
    j=rank;
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    j=j2;
    MPI_Bcast(&begblock[j2],1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  j2=begblock[0];
  begblock[0]=0;
  for(i=1; i<mpisize; i++) {
    j1=begblock[i];
    begblock[i]=begblock[i-1]+j2;
    j2=j1;
  }
  //for(i=0; i<mpisize; i++)printf("rank %d beg block %d\n",rank,begblock[i]);
  j2=nmatin;
  MPI_Bcast(&j2,1, MPI_INT,0, PETSC_COMM_WORLD);
  ndblockinrank[0]=j2;
  for(i=1; i<mpisize; i++) {
    j2=nmatin;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    ndblockinrank[i]=ndblockinrank[i-1]+j2;
  }
  //int *birank= (int *) calloc (nmatin,sizeof(int));
  //int *cirank= (int *) calloc (nmatin,sizeof(int));
  //for(i=0; i<nmatin; i++)cirank[i]=-1;
  //int *aord= (int *) calloc (nmatin,sizeof(int));
  for(i=0; i<ndblock; i++) {
    j2=i-ndblockinrank[rank];
    if(j2<0&&j2>=-nmatin) {
      j1=i/(nreg+1);
      j=i-(nreg+1)*j1;
      if(rank==0)j3=0;
      else j3=ndblockinrank[rank-1];
      //if(j<nreg)cirank[j2+nmatin]=(nreg+1)*j1+nreg-j3;
      //else cirank[j2+nmatin]=-1;
    }
  }
  j1=0;
  for(j=0; j<mpisize; j++)if(j<rank)j1+=ndblockinrank[j];
  free(ndblockinrank);
//   for(i=0; i<nmatin; i++) {
//     j2=(j1+i)/(nreg+1);
//     if(i-j2*(nreg+1)!=0)aord[i]=1;
//   }
  ierr = PetscMalloc(nmatin*sizeof(IS **),&rowindices);
  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*sizeof(IS **),&rowindicesbc);
//  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&colindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbc1);
  CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbc2);
  //CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*sizeof(IS **),&Cindices);
//  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*sizeof(IS **),&Bindices);
//  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatinBB*sizeof(IS **),&BBindices);
//  CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatAij);
  //CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatint*sizeof(Mat *),&submatC);
  //CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatinBB*sizeof(Mat *),&submatBB);
//  CHKERRQ(ierr);
  Mat** submatBij = (Mat**)calloc(nmatin,sizeof(Mat*));
//   for (i=0; i<nmatin; i++) {
//     submatBij[i] = (Mat*)calloc(1,sizeof(Mat));
//   }
  Mat** submatCij = (Mat**)calloc(nmatin,sizeof(Mat*));
//   for (i=0; i<nmatin; i++) {
//     submatCij[i] = (Mat*)calloc(1,sizeof(Mat));
//   }
//  Mat** submatB1 = (Mat**)calloc(nmatint,sizeof(Mat*));
//  for (i=0; i<nmatint; i++) {
//    submatB1[i] = (Mat*)calloc(1,sizeof(Mat));
//  }
//  Mat** submatB2 = (Mat**)calloc(nmatint,sizeof(Mat*));
//  for (i=0; i<nmatint; i++) {
//    submatB2[i] = (Mat*)calloc(1,sizeof(Mat));
//  }

//  ha_cgetype **yi2= (ha_cgetype**)calloc(nmatin,sizeof(ha_cgetype*));
//  for (i=0; i<nmatin; i++) {
//    yi2[i] = (ha_cgetype*)calloc(ha_ndblocks[i+begblock[rank]],sizeof(ha_cgetype));
//  }
//  sumrowcolin=0;
//  for(j=0; j<ndblock; j++)sumrowcolin+=ha_ndblocks[j];
//  PetscInt *offblock= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
//  PetscInt *offblockrow= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
//  offblock[0]=0;
//  for(j=0; j<ndblock; j++)offblock[j+1]=offblock[j]+countvarintra1[j+1]-countvarintra1[j]-ha_ndblocks[j];
//  offblockrow[0]=0;
//  for(j=0; j<ndblock; j++)offblockrow[j+1]=offblockrow[j]+counteqnoadd[j]-ha_ndblocks[j];

  char **fn01= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn01[i] = (char*)calloc(1024,sizeof(char));
  char **fn02= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn02[i] = (char*)calloc(1024,sizeof(char));
  char **fn03= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn03[i] = (char*)calloc(1024,sizeof(char));

  //timestr=clock();
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      if(j1<10)strcpy(j1name,"000");
      if(j1<100&&j1>9)strcpy(j1name,"00");
      if(j1<1000&&j1>99)strcpy(j1name,"0");
      if(j1>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j1);
      strcat(j1name,filename);
      strcpy(filename,"_irnv");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn01[j1],filename);
      strcpy(filename,"_keep");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn02[j1],filename);
      strcpy(filename,"_vav");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn03[j1],filename);
    }
  }
  for(j1=0; j1<nmatinplus; j1++) {
    free(fn01[j1]);//1
    free(fn02[j1]);//1
    free(fn03[j1]);//1
  }
  free(fn01);//1
  free(fn02);//1
  free(fn03);//1
  //printf("OK012345 rank %d mpis %d sizeof int %ld size of long int %ld\n",rank,mpisize,sizeof(int),sizeof(long int));
//  for(j=0; j<mpisize; j++) {
//    //printf("j %d\n",j);
//    for(j3=0; j3<nmatinplus; j3++) {
//      if(j3>=nmatin)i=nmatin-1;
//      else i=j3;
//      //printf("rank %d OK1!!! i %d j %d nmatin %d\n",rank,i,j,nmatin);
//      bfirst=counteq[i+begblock[j]];
//      bend=ha_ndblocks[i+begblock[j]];
//      MPI_Bcast(&bfirst,1, MPI_INT,j, PETSC_COMM_WORLD);
//      MPI_Bcast(&bend,1, MPI_INT,j, PETSC_COMM_WORLD);
//      PetscInt *indices= (PetscInt *) calloc (bend,sizeof(PetscInt));
//      ha_cgetype *yi0 = (ha_cgetype*)calloc(bend,sizeof(ha_cgetype));
//      for(j1=0; j1<bend; j1++) { //ha_ndblocks[i+begblock[j]]
//        j2=ha_rows[bfirst+j1];//+bfirst;
//        if(j2>=Istart&&j2<Iend)indices[j1]=j2;
//        else indices[j1]=-1;
//      }
//      //printf("rank %d i %d j %d bf %d be %d\n",rank,i,j,bfirst,bend);
//      VecGetValues(b,bend,indices,yi0);
//      HaReduce(yi0,(forint)bend,mpisize,rank,j);
//      if(rank==j&&i<nmatin) {
//        memcpy (yi2[i],yi0,bend*sizeof(ha_cgetype));
//      }
//      free(indices);//1
//      free(yi0);//1
//    }
//  }
//  ha_cgetype **yi1= (ha_cgetype**)calloc(nmatint,sizeof(ha_cgetype*));
//  for (i=0; i<nmatint; i++) {
//    j1=0;
//    for (j=0; j<nreg+1; j++) {
//      j1+=ha_ndblocks[i*(nreg+1)+j+begblock[rank]];
//    }
//    yi1[i] = (ha_cgetype*)calloc(j1,sizeof(ha_cgetype));
//  }
//  for (i=0; i<nmatint; i++) {
//    j1=0;
//    for (j=0; j<nreg+1; j++) {
//      memcpy (&yi1[i][j1],yi2[i*(nreg+1)+j],ha_ndblocks[i*(nreg+1)+j+begblock[rank]]*sizeof(ha_cgetype));
//      free(yi2[i*(nreg+1)+j]);//1
//      j1+=ha_ndblocks[i*(nreg+1)+j+begblock[rank]];
//    }
//  }
//  free(yi2);//1
  //printf("OK12345 rank %d\n",rank);
  j1=0;
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (ha_ndblocks[i+begblock[rank]],sizeof(PetscInt));
    //printf("ha_row 0 %d 588 %d\n",ha_rows[0],ha_rows[588]);
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_rows[bfirst+j];//+bfirst;
    }

    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
    bfirst=countvarintra1[i+begblock[rank]];
    //if(rank==0)printf("rank %d i %d bfirst %d block %d\n",rank,i,bfirst,ha_ndblocks[i+begblock[rank]]);
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_cols[bfirst+j];//+bfirst;
    }
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,colindices+i);
    free(indices);//1
  }
  for(i=0; i<nmatint; i++) {
    j2=0;
    for(j1=0; j1<nreg+1; j1++)j2+=ha_ndblocks[i*(nreg+1)+j1+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (j2,sizeof(PetscInt));
//    j2=0;
//    for(j1=0; j1<nreg+1; j1++) {
//      bfirst=counteq[i*(nreg+1)+j1+begblock[rank]];
//      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
//        indices[j2]=ha_rows[bfirst+j];//+bfirst;
//        j2++;
//      }
//    }
//    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,rowindicesbc+i);
    j2=0;
    for(j1=0; j1<nreg; j1++) {
      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
        indices[j2]=ha_cols[bfirst+j];//+bfirst;
        j2++;
      }
    }
    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbc1+i);
//    j2=0;
//    for(j1=nreg; j1<nreg+1; j1++) {
//      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
//      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
//        indices[j2]=ha_cols[bfirst+j];//+bfirst;
//        j2++;
//      }
//    }
    //ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbc2+i);
    free(indices);
  }
//  printf("0okkkk\n");
//  PetscInt *indicesC= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
//  PetscInt *indicesB= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
//  for(i=0; i<ndblock; i++) {
//    bend=countvarintra1[i+1];
//    for(j=offblock[i]; j<offblock[i+1]; j++) {
//      j1=j-offblock[i];
//      indicesC[j]=ha_cols[bend-j1-1];//+bfirst;
//    }
//  }
//  j1=offblock[ndblock];
//  for(i=countvarintra1[ndblock]; i<VecSize; i++) {
//    indicesC[j1]=i;
//    j1++;
//  }
//  printf("1okkkk\n");
//  for(i=0; i<ndblock; i++) {
//    bend=counteq[i]+counteqnoadd[i];
//    for(j=offblockrow[i]; j<offblockrow[i+1]; j++) {
//      j1=j-offblockrow[i];
//      indicesB[j]=ha_rows[bend-j1-1];//+bfirst;
//    }
//  }
//  j1=offblockrow[ndblock];
//  for(i=counteq[ndblock-1]+counteqnoadd[ndblock-1]; i<VecSize; i++) {
//    indicesB[j1]=i;
//    j1++;
//  }
//  printf("2okkkk\n");
//  printf("rank1 %d\n",rank);
//  for(i=0; i<VecSize-sumrowcolin; i++)if(indicesB[i]>=Istart) {
//      j=i;
//      break;
//    }
//  for(j2=i; j2<VecSize-sumrowcolin; j2++)if(indicesB[j2]>=Iend) {
//      j1=j2;
//      break;
//    }
//  if(j2==VecSize-sumrowcolin)j1=VecSize-sumrowcolin;
//  printf("rank %d j %d j1 %d istart %d iend %d\n",rank,j,j1,Istart,Iend);
//  ha_cgetype *vecbiui= (ha_cgetype *) calloc (VecSize-sumrowcolin,sizeof(ha_cgetype));
//  VecGetValues(b,j1-j,&indicesB[j],&vecbiui[j]);//implicite yd save mem
  //ierr = VecDestroy(&b);
  //CHKERRQ(ierr);
//  printf("rank2 %d\n",rank);
//  for(i=0; i<nmatint; i++) {
//    ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesC,PETSC_COPY_VALUES,Cindices+i);
//    ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesB,PETSC_COPY_VALUES,Bindices+i);
//  }
//  free(indicesB);//1
//  free(offblock);//1
//  free(offblockrow);//1
//  PetscInt *indexBB= (PetscInt *) calloc (VecSize,sizeof(PetscInt));
//  for(i=0; i<VecSize; i++)indexBB[i]=i;
//  ISCreateGeneral(PETSC_COMM_SELF,VecSize,indexBB,PETSC_COPY_VALUES,BBindices);
  printf("Begin preparation rank %d\n",rank);
  ierr = MatCreateSubMatrices(A,nmatin,rowindices,colindices,MAT_INITIAL_MATRIX,&submatAij);
  //if(rank==0)printf("rank %d i %d bbrowij %d\n",rank,(nmatint-1)*(nreg+1)+nreg,submatAij[(nmatint-1)*(nreg+1)+nreg]->cmap->n);
  int insizes=17;
  int *insize=(int *) calloc (insizes*nmatin,sizeof(int));
  printf("rank %d A\n",rank);
  PetscInt *indicesbbij=NULL;//(int *) calloc (1,sizeof(PetscInt));
  ierr = PetscMalloc(1*sizeof(IS **),&rowBBij);
  CHKERRQ(ierr);
  ierr = PetscMalloc(1*sizeof(IS **),&colBBij);
  CHKERRQ(ierr);

  for(i=0; i<nmatminust; i++) {
    j2=i*(nreg+1)+nreg;
    ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindicesbc1+i,MAT_INITIAL_MATRIX,&submatBBij);
    nrow=submatBBij[0]->rmap->n;
    ncol=submatBBij[0]->cmap->n;
    if(nrow>ncol) {
      indicesbbij=realloc(indicesbbij,nrow*sizeof(PetscInt));
      for(j=0; j<nrow; j++)indicesbbij[j]=j;
    }
    else {
      indicesbbij=realloc(indicesbbij,ncol*sizeof(PetscInt));
      for(j=0; j<ncol; j++)indicesbbij[j]=j;
    }
    ISCreateGeneral(PETSC_COMM_SELF,nrow,indicesbbij,PETSC_COPY_VALUES,rowBBij);
    j3=0;
    for(j=0; j<nreg; j++) {
      j1=i*(nreg+1)+j;
      ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
      ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[j1+begblock[rank]],indicesbbij+j3,PETSC_COPY_VALUES,colBBij);
      ierr = MatCreateSubMatrices(submatBBij[0],1,rowBBij,colBBij,MAT_INITIAL_MATRIX,&submatBij[j1]);
      ierr = ISDestroy(&colBBij[0]);
      j3+=ha_ndblocks[j1+begblock[rank]];
      //printf("rank %d i %d j3 %d bbrowij %d\n",rank,i,j,submatCij[j1][0]->cmap->n);
    }
    //printf("i %d\n",i);
    ierr = MatDestroy(&submatBBij[0]);
    ierr = PetscFree(submatBBij);//1
    submatBBij=NULL;
    ierr = ISDestroy(&rowBBij[0]);
  }

//  for(i=0; i<nmatminust; i++) {
//      j2=i*(nreg+1)+nreg;
//    for(j=0; j<nreg; j++) {
//      j1=i*(nreg+1)+j;
////      ierr = ISSort(rowindices[j2]);
////      ierr = ISSort(rowindices[j1]);
////      ierr = ISSort(colindices[j1]);
////      ierr = ISSort(colindices[j2]);
//      ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
//      //printf("i %d j %d C\n",i,j);
//      ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindices+j1,MAT_INITIAL_MATRIX,&submatBij[j1]);
//      //printf("i %d j %d B\n",i,j);
//      //printf("rank %d i %d j3 %d bbrowij %d\n",rank,i,j,submatCij[j1][0]->cmap->n);
//    }
//  }
//  //printf("rank %d A1\n",rank);
  if(unequal==1){
  if(nmatint<nmatinplust) for(i=0; i<nreg; i++) {//to create again! MatCreateSubMatrices must be the same across
      MatDestroy(&submatCij[(nmatint-1)*(nreg+1)+i][0]);
      ierr = PetscFree(submatCij[(nmatint-1)*(nreg+1)+i]);
      MatDestroy(&submatBij[(nmatint-1)*(nreg+1)+i][0]);
      ierr = PetscFree(submatBij[(nmatint-1)*(nreg+1)+i]);
    }
  j2=(nmatint-1)*(nreg+1)+nreg;
  ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindicesbc1+(nmatint-1),MAT_INITIAL_MATRIX,&submatBBij);
  nrow=submatBBij[0]->rmap->n;
  ncol=submatBBij[0]->cmap->n;
  if(nrow>ncol) {
    indicesbbij=realloc(indicesbbij,nrow*sizeof(PetscInt));
    for(j=0; j<nrow; j++)indicesbbij[j]=j;
  }
  else {
    indicesbbij=realloc(indicesbbij,ncol*sizeof(PetscInt));
    for(j=0; j<ncol; j++)indicesbbij[j]=j;
  }
  ISCreateGeneral(PETSC_COMM_SELF,nrow,indicesbbij,PETSC_COPY_VALUES,rowBBij);
  j3=0;
  for(j=0; j<nreg; j++) {
    j1=(nmatint-1)*(nreg+1)+j;
    ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[j1+begblock[rank]],indicesbbij+j3,PETSC_COPY_VALUES,colBBij);
    ierr = MatCreateSubMatrices(submatBBij[0],1,rowBBij,colBBij,MAT_INITIAL_MATRIX,&submatBij[j1]);
    ierr = ISDestroy(&colBBij[0]);
    j3+=ha_ndblocks[j1+begblock[rank]];
    //printf("rank %d i %d j3 %d bbrowij %d\n",rank,i,j,submatCij[j1][0]->cmap->n);
  }
  ierr = MatDestroy(&submatBBij[0]);
  ierr = PetscFree(submatBBij);//1
  submatBBij=NULL;
  ierr = ISDestroy(&rowBBij[0]);
  }
  free(begblock);
  free(indicesbbij);
  PetscFree(rowBBij);
  PetscFree(colBBij);

//  for(j=0; j<nreg; j++) {
//    j1=(nmatint-1)*(nreg+1)+j;
//    j2=(nmatint-1)*(nreg+1)+nreg;
//    ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
//    ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindices+j1,MAT_INITIAL_MATRIX,&submatBij[j1]);
//    //printf("rank %d i %d j3 %d bbrowij %d\n",rank,i,j1,submatCij[j1][0]->cmap->n);
//  }
  //printf("rank %d A\n",rank);
  //ierr = MatCreateSubMatrices(A,nmatint,rowindicesbc,Cindices,MAT_INITIAL_MATRIX,&submatC);
  //printf("C\n");
//  ierr = MatCreateSubMatrices(A,nmatinBB,Bindices,BBindices,MAT_INITIAL_MATRIX,&submatBB);
  //ierr = MatDestroy(&A);
  //CHKERRQ(ierr);
//  PetscInt BBrow,BBcol;//,CCrow,CCcol;
//  MatGetSize(submatBB[0],&BBrow,&BBcol);
//  //MatGetSize(submatC[0],&CCrow,&CCcol);
//  printf("rank %d bbrow %d bbcol %d\n",rank,BBrow,BBcol);//,CCrow,CCcol);
//  PetscInt *indexBBi= (PetscInt *) calloc (BBrow,sizeof(PetscInt));
//  for(i=0; i<BBrow; i++)indexBBi[i]=i;
//  for(i=0; i<nmatint; i++) {
//    ierr = ISDestroy(&Bindices[i]);
//    ISCreateGeneral(PETSC_COMM_SELF,BBrow,indexBBi,PETSC_COPY_VALUES,Bindices+i);
//  }
  //ierr = MatCreateSubMatrices(submatBB[0],nmatint,Bindices,colindicesbc,MAT_INITIAL_MATRIX,&submatB);
  //char           ordering[256] = MATORDERINGRCM;
//  for(i=0; i<nmatint; i++) {
//    ierr = MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc1[i],MAT_INITIAL_MATRIX,&submatB1[i][0]);
//    ierr = ISDestroy(&colindicesbc1[i]);
//    ierr = MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc2[i],MAT_INITIAL_MATRIX,&submatB2[i][0]);
//    ierr = ISDestroy(&colindicesbc2[i]);
//  }
//  ierr = PetscFree(colindicesbc1);//1
//  ierr = PetscFree(colindicesbc2);//1
//  printf("B\n");
//  if(rank==mpisize-1)proc1=1;
//  ierr = PetscMalloc(proc1*sizeof(Mat *),&submatD);
//  CHKERRQ(ierr);
//  printf("D\n");
//  MatCreateSubMatrices(submatBB[0],proc1,Bindices,Cindices,MAT_INITIAL_MATRIX,&submatD);
  //printf("End partitionning rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&rowindices[i]);//1
    CHKERRQ(ierr);
    ierr = ISDestroy(&colindices[i]);//1
    CHKERRQ(ierr);
  }
  for (i=0; i<nmatint; i++) {
    ierr = ISDestroy(&colindicesbc1[i]);//1
    CHKERRQ(ierr);
  }
//  for (i=0; i<nmatint; i++) {
//    ierr = ISDestroy(&Cindices[i]);//1
//    CHKERRQ(ierr);
//    ierr = ISDestroy(&Bindices[i]);//1
//    CHKERRQ(ierr);
//  }
  PetscFree(colindices);//1
//  ierr = ISDestroy(&BBindices[0]);//1
  PetscFree(rowindices);//1
  PetscFree(colindicesbc1);//1

//  for (i=0; i<nmatint; i++) {
//    ierr = ISDestroy(&rowindicesbc[i]);//1
//    CHKERRQ(ierr);
//  }
//  PetscFree(rowindicesbc);//1

//  PetscFree(Cindices);//1
//  PetscFree(Bindices);//1
//  PetscFree(BBindices);//1
//  PetscFree(indexBB);//1
  //PetscFree(indexBBi);//1
//  ierr = MatDestroy(&submatBB[0]);//1
//  PetscFree(submatBB);//1
//  CHKERRQ(ierr);
//  printf("OK1 rank %d sumrowcolin %d\n",rank,sumrowcolin);
  //ha_cgetype *xi1 = (ha_cgetype*)calloc(sumrowcolin,sizeof(ha_cgetype));
  //ccolnorms=realloc(ccolnorms,1*sizeof(PetscReal));
  int bivirowsize=1,bivicolsize=1,bbrowij,ddrowi;
  long int vecbivisize=0,li,lj;//,halfvec;
  //ha_cgetype *xi1point;
  //ha_cgetype *vecbivi= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
  //uvadd xi1indx=0;
  long int nz0,nz1,nz2,nz3=j,halfj2;
  long int lj2;
  int jthrd,nthrd;
  //omp_set_num_threads(1);//shut off openmp here due to error
  //FILE** npresolfilea=(FILE**) calloc (nmatint,sizeof(FILE*));
  //FILE** npresolfilea1=(FILE**) calloc (nmatint,sizeof(FILE*));
  //FILE** npresolfileb=(FILE**) calloc (nmatint,sizeof(FILE*));
  //FILE** npresolfileb1=(FILE**) calloc (nmatint,sizeof(FILE*));
  //FILE** npresolfilec=(FILE**) calloc (nmatint,sizeof(FILE*));
  //FILE** npresolfilec1=(FILE**) calloc (nmatint,sizeof(FILE*));
  //omp_set_num_threads(1);
  #pragma omp parallel private(jthrd,nthrd,j3,j4,bivirowsize,bivicolsize,bbrowij,ai,nz,nrow,i,j,j1,j2,li,lj,ddrowi,vecbivisize,aic,ajc,valsc,nzc,nrowc,ncolc,ncolb,nrowb,aj,vals,lj2,nz0,j1name,filename,presolfile,fwrt,fd1,nz1,cntl6in,ncol,lasize,ldsize) shared(insize,submatAij,submatBij,submatCij)
  //,npresolfilea,npresolfilea1,npresolfileb,npresolfileb1,npresolfilec,npresolfilec1,timestr
  {
  long int *bivinzrow=NULL,irnmems=0;//(long int *) calloc (1,sizeof(long int));
  PetscInt *bivinzcol=NULL;//(PetscInt *) calloc (1,sizeof(PetscInt));
  PetscReal *ccolnorms=NULL;//(PetscReal *) calloc (1,sizeof(PetscReal));
  long int *bivinzrow1=NULL;//(long int *) calloc (1,sizeof(long int));
  long int *bivinzcol1=NULL;//(long int *) calloc (1,sizeof(long int));
  int *irn=NULL,*jcn=NULL,*irn1=NULL,*jcn1=NULL;
  int windx=0,bindx,eindx;
  //size_t *nfwrt=(size_t *) calloc (nmatint,sizeof(size_t));
  jthrd=omp_get_thread_num();
  nthrd=omp_get_max_threads( );
  nthrd=90/nthrd;
  bindx=nthrd*jthrd;
  eindx=bindx+nthrd;
  windx=bindx;
  Mat_SeqAIJ         *aa,*ac,*ab;
  Mat submatCT;
  ha_cgetype *vecbivi=NULL;
  long int bbrowijmax=0,bivirowsizemax=0,bivicolsizemax=0,vecbivisizemax=0,lasizemax=0,nz1max=0,nrowmax=0,ncolmax=0,nrowamax=0,ncolamax=0,nebmax=0;
  #pragma omp for schedule (static)
  for(j3=0; j3<nmatint; j3++) {
    bbrowij=submatBij[j3*(nreg+1)][0]->rmap->n;
    if(bbrowij>bbrowijmax)bbrowijmax=bbrowij;
  }
  ccolnorms=realloc(ccolnorms,bbrowijmax*sizeof(PetscReal));
  bivinzrow=realloc(bivinzrow,bbrowijmax*sizeof(long int));
  bivinzcol=realloc(bivinzcol,bbrowijmax*sizeof(PetscInt));
  #pragma omp for schedule (static)
  for(j3=0; j3<nmatint; j3++) {
    bivirowsize=1;
    bivicolsize=1;
    bbrowij=submatBij[j3*(nreg+1)][0]->rmap->n;
    memset(bivinzrow,0,bbrowij*sizeof(long int));
    memset(bivinzcol,0,bbrowij*sizeof(PetscInt));
    for(j1=0; j1<nreg; j1++) {
      MatGetColumnNorms(submatCij[j1+j3*(nreg+1)][0],NORM_1,ccolnorms);
      for(i=0; i<bbrowij; i++)if(ccolnorms[i]>0)bivinzcol[i]++;
      ab=(Mat_SeqAIJ*)submatBij[j1+j3*(nreg+1)][0]->data;//*aa=subA->data;
      ai= ab->i;
      nz=ab->nz;
      if(nebmax<nz)nebmax=nz;
      nrow=submatBij[j1+j3*(nreg+1)][0]->rmap->n;
      //printf("j1 %d nrow %d\n",j1,nrow);
      for(i=0; i<nrow-1; i++)if(ai[i]!=ai[i+1])bivinzrow[i]++;
      if(ai[i]<nz)bivinzrow[i]++;
    }
    ddrowi=submatBij[j3*(nreg+1)][0]->rmap->n;
    for(i=0; i<ddrowi; i++) {
      if(bivinzrow[i]>0) {
        bivinzrow[i]=bivirowsize;
        bivirowsize++;
      }
      if(bivinzcol[i]>0) {
        bivinzcol[i]=bivicolsize;
        bivicolsize++;
      }
    }
    bivirowsize--;
    bivicolsize--;
    if(bivirowsize>bivirowsizemax)bivirowsizemax=bivirowsize;
    if(bivicolsize>bivicolsizemax)bivicolsizemax=bivicolsize;
    vecbivisize=((long int)bivirowsize)*((long int)bivicolsize);
    if(vecbivisize>vecbivisizemax)vecbivisizemax=vecbivisize;
    for(j1=0; j1<nreg; j1++) {
      j4=j1+j3*(nreg+1);
      aa=(Mat_SeqAIJ*)submatAij[j4]->data;//*aa=subA->data;
      nrow=submatAij[j4]->rmap->n;
      if(nrowamax<nrow)nrowamax=nrow;
      ncol=submatAij[j4]->cmap->n;
      if(ncolamax<ncol)ncolamax=ncol;
      nz=aa->nz;
      lasize=ceil((laA/100.0)*nz);
      if(lasize>lasizemax)lasizemax=lasize;
    }
  }
  lasizemax+=10;

      int *insized=(int *) calloc (5+nreg*insizes,sizeof(int));
    bivinzrow1=realloc (bivinzrow1,bivirowsizemax*sizeof(long int));
    if(bivinzrow1==NULL)printf("Error!!!! No memory left for bivinzrow1!\n");
    bivinzcol1=realloc (bivinzcol1,bivicolsizemax*sizeof(long int));
    if(bivinzcol1==NULL)printf("Error!!!! No memory left for bivinzcol1!\n");
      vecbivi= (ha_cgetype *) malloc (vecbivisizemax*sizeof(ha_cgetype));
    if(vecbivi==NULL)printf("Error!!!! No memory left for vecbivi!\n");
      irn=(int *) malloc (lasizemax*sizeof(int));
    if(irn==NULL)printf("Error!!!! No memory left for irn!\n");
      jcn=(int *) malloc (lasizemax*sizeof(int));
    if(jcn==NULL)printf("Error!!!! No memory left for jcn!\n");
      ha_cgetype *values= (ha_cgetype *) calloc (lasizemax,sizeof(ha_cgetype));
    if(values==NULL)printf("Error!!!! No memory left for values!\n");
      int* jcnb1=(int *) malloc (nebmax*sizeof(int));
    if(jcnb1==NULL)printf("Error!!!! No memory left for jcnb1!\n");
      ha_cgetype *sol48= (ha_cgetype *) malloc (nrowamax*sizeof(ha_cgetype));
    if(sol48==NULL)printf("Error!!!! No memory left for sol48!\n");
      ha_cgetype *b48= (ha_cgetype *) malloc (nrowamax*sizeof(ha_cgetype));
    if(b48==NULL)printf("Error!!!! No memory left for b48!\n");
      int* keep=NULL,*iw51=NULL;
  keep=(int *) malloc ((nrowamax+5*ncolamax+4*ncolamax+7)*sizeof(int));//ICNTL(6)=1
    if(keep==NULL)printf("Error!!!! No memory left for keep!\n");
  ha_cgetype *w51=NULL;
  if(nrowamax>ncolamax) w51=(ha_cgetype *) malloc ((5*nrowamax)*sizeof(ha_cgetype));
  else w51=(ha_cgetype *) malloc ((5*ncolamax)*sizeof(ha_cgetype));
    if(w51==NULL)printf("Error!!!! No memory left for w51!\n");
  iw51=(int *) malloc ((6*nrowamax+3*ncolamax)*sizeof(int));
    if(iw51==NULL)printf("Error!!!! No memory left for iw51!\n");
  #pragma omp for schedule (static)
  for(j3=0; j3<nmatint; j3++) {
    bivirowsize=1;
    bivicolsize=1;
    bbrowij=submatBij[j3*(nreg+1)][0]->rmap->n;
    printf("j3 %d bbrowij %djth %d windx %d\n",j3,bbrowij,jthrd,windx);
    //ccolnorms=realloc(ccolnorms,bbrowij*sizeof(PetscReal));
    //bivinzrow=realloc(bivinzrow,bbrowij*sizeof(long int));
    memset(bivinzrow,0,bbrowij*sizeof(long int));
    //bivinzcol=realloc(bivinzcol,bbrowij*sizeof(PetscInt));
    memset(bivinzcol,0,bbrowij*sizeof(PetscInt));
    for(j1=0; j1<nreg; j1++) {
      MatGetColumnNorms(submatCij[j1+j3*(nreg+1)][0],NORM_1,ccolnorms);
      for(i=0; i<bbrowij; i++)if(ccolnorms[i]>0)bivinzcol[i]++;
      ab=(Mat_SeqAIJ*)submatBij[j1+j3*(nreg+1)][0]->data;//*aa=subA->data;
      ai= ab->i;
      nz=ab->nz;
      nrow=submatBij[j1+j3*(nreg+1)][0]->rmap->n;
      //printf("j1 %d nrow %d\n",j1,nrow);
      for(i=0; i<nrow-1; i++)if(ai[i]!=ai[i+1])bivinzrow[i]++;
      if(ai[i]<nz)bivinzrow[i]++;
    }
    ddrowi=submatBij[j3*(nreg+1)][0]->rmap->n;
    //printf("ddrowi %d\n",ddrowi);
    for(i=0; i<ddrowi; i++) {
      if(bivinzrow[i]>0) {
        bivinzrow[i]=bivirowsize;
        bivirowsize++;
      }
      if(bivinzcol[i]>0) {
        bivinzcol[i]=bivicolsize;
        bivicolsize++;
      }
    }
    //printf("j3 %d\n",j3);
    bivirowsize--;
    bivicolsize--;
    //bivinzrow1=realloc (bivinzrow1,bivirowsize*sizeof(long int));
    //bivinzcol1=realloc (bivinzcol1,bivicolsize*sizeof(long int));
    j=0;
    j2=0;
    for(li=0; li<ddrowi; li++) {
      if(bivinzrow[li]>0) {
        bivinzrow1[j]=li*bbrowij;
        j++;
      }
      if(bivinzcol[li]>0) {
        bivinzcol1[j2]=li;
        j2++;
      }
      bivinzrow[li]=(bivinzrow[li]-1)*bivicolsize;//faster MSOL
    }
    vecbivisize=((long int)bivirowsize)*((long int)bivicolsize);
    printf("rank %d vecbivisize %ld bivirowsize %d bivicolsize %d\n",rank,vecbivisize,bivirowsize,bivicolsize);
    if(vecbivisize<0)printf("Overflow detected, the interface problem is too big!!!!!!\n");
    //vecbivi=realloc (vecbivi,vecbivisize*sizeof(int));
    //vecbivi= (ha_cgetype *) calloc (vecbivisize,sizeof(ha_cgetype));
    memset(vecbivi,0,vecbivisize*sizeof(ha_cgetype));
    //printf("vecbivisize %ld rank %d\n",vecbivisize,rank);
    //xi1indx=0;
    //printf("rank0 %d j4 %d nmatin %d\n",rank,j4,nmatin);
//       irn=(int *) malloc (lasizemax*sizeof(int));
//       jcn=(int *) malloc (lasizemax*sizeof(int));
//       ha_cgetype *values= (ha_cgetype *) malloc (lasizemax*sizeof(ha_cgetype));
//       int* jcnb1=(int *) malloc (nebmax*sizeof(int));
//       ha_cgetype *sol48= (ha_cgetype *) malloc (nrowamax*sizeof(ha_cgetype));
//       ha_cgetype *b48= (ha_cgetype *) malloc (nrowamax*sizeof(ha_cgetype));
//   int *keep=(int *) malloc ((nrowamax+5*ncolamax+4*ncolamax+7)*sizeof(int));//ICNTL(6)=1
//   ha_cgetype *w51;
//   if(nrowamax>ncolamax) w51=(ha_cgetype *) malloc ((5*nrowamax)*sizeof(ha_cgetype));
//   else w51=(ha_cgetype *) malloc ((5*ncolamax)*sizeof(ha_cgetype));
//   int *iw51=(int *) malloc ((6*nrowamax+3*ncolamax)*sizeof(int));
    //#pragma omp for
    for(j1=0; j1<nreg; j1++) {
      j4=j1+j3*(nreg+1);
      //time(&timestr);//=clock();
      //printf("rank %d j4 %d nmatin %d\n",rank,j4,nmatin);
      aa=(Mat_SeqAIJ*)submatAij[j4]->data;//*aa=subA->data;
      MatTranspose(submatCij[j4][0],MAT_INITIAL_MATRIX,&submatCT);
      MatDestroy(&submatCij[j4][0]);//1
      ac=(Mat_SeqAIJ*)submatCT->data;//*aa=subA->data;
      ai= aa->i;
      aj= aa->j;
      vals=aa->a;
      nz=aa->nz;
      nrow=submatAij[j4]->rmap->n;
      ncol=submatAij[j4]->cmap->n;
      lasize=ceil((laA/100.0)*nz);
      //printf("la %ld lam %ld\n",lasize,lasizemax);
//       irn=(int *) calloc (lasize,sizeof(int));
//       jcn=(int *) calloc (lasize,sizeof(int));
//       values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
      //memset(irn,0,lasize*sizeof(int));
      //memset(jcn,0,lasize*sizeof(int));
      //memset(values,0,lasize*sizeof(ha_cgetype));
      insize[j4*insizes+13]=bivirowsize;
      insize[j4*insizes+14]=bivicolsize;
      insize[j4*insizes]=nrow;
      insize[j4*insizes+1]=ncol;
      insize[j4*insizes+2]=nz;
      insize[j4*insizes+15]=7+windx;//7+rank+jthrd+windx;//(j4+jthrd)%90+7;
      windx++;
      if(windx==eindx)windx=bindx;
      aic= ac->i;
      ajc= ac->j;
      valsc=ac->a;
      nzc=ac->nz;
      nrowc=submatCT->rmap->n;
      ncolc=submatCT->cmap->n;
      insize[j4*insizes+3]=nrowc;
      insize[j4*insizes+4]=ncolc;
      memcpy (irn+nz,ai,(nrow+1)*sizeof(PetscInt));
      memcpy (jcn,aj,nz*sizeof(PetscInt));
      memcpy (values,vals,nz*sizeof(ha_cgetype));
      //printf("nz %d j %d\n",nz,j);
      insize[j4*insizes+2]=nz;
      insize[j4*insizes+5]=nzc;
      MatDestroy(&submatAij[j4]);//1
      insize[j4*insizes+9]=laA;
      insize[j4*insizes+16]=lasize;
      insize[j4*insizes+10]=rank;
      insize[j4*insizes+11]=j4;
      MatGetSize(submatBij[j4][0],&nrowb,&ncolb);
      ab=(Mat_SeqAIJ*)submatBij[j4][0]->data;//*aa=subA->data;
      ai= ab->i;//irnbs+j1irnbs;
      aj= ab->j;//jcnbs+j1nz;
      vals=ab->a;//valbs+j1nz;
      nz=ab->nz;//nzbs[j1];
      insize[j4*insizes+6]=nrowb;//ncolb;
      insize[j4*insizes+7]=ncolb;//nrowb;
      insize[j4*insizes+8]=nz;
      prep48m_msol_(insize+j4*insizes,irn,jcn,values,aic,ajc,valsc,ai,aj,vals,vecbivi,bivinzrow,bivinzcol,jcnb1,sol48,b48,w51,iw51,keep);
      insize[j4*insizes+15]=0;
      MatDestroy(&submatCT);
      MatDestroy(&submatBij[j4][0]);//1
      //free(yi1[j1]);
      //time(&timeend);
      //printf("t %d Submatrix %d rank %d calculation time %f\n",j3,j1,rank,difftime(timeend,timestr));
    }
//       free(irn);//1
//       irn=NULL;
//       free(jcn);//1
//       jcn=NULL;
//       free(values);//1
//     free(jcnb1);
//     jcnb1=NULL;
//     free(sol48);
//     sol48=NULL;
//     free(b48);
//     b48=NULL;
//     free(keep);
//     keep=NULL;
//     free(w51);
//     w51=NULL;
//     free(iw51);
//     iw51=NULL;

    j4=nreg+j3*(nreg+1);
    aa=(Mat_SeqAIJ*)submatAij[j4]->data;//*aa=subA->data;
    ai= aa->i;
    aj= aa->j;
    vals=aa->a;
    nz=aa->nz;

    lj=0;
    for(li=vecbivisize-1; li!=-1; --li)if(vecbivi[li]!=0)lj++;
    irnmems=nz+lj;
//     ldsize=ceil((laDi/100.0)*irnmems);//first estimate
//     if(vecbivisize<ldsize){
//       vecbivi=realloc(vecbivi,ldsize*sizeof(ha_cgetype));
//       //vecbivisize=ldsize;
//     }
//     if(ldsize>irnmems)irnmems=ldsize;
    #pragma omp critical
    {
    irn1=realloc(irn1,irnmems*sizeof(int));
    if(irn1==NULL)printf("Error!!!! No memory left for irn1!\n");
    jcn1=realloc(jcn1,irnmems*sizeof(int));
    if(jcn1==NULL)printf("Error!!!! No memory left for jcn1!\n");
//     int *irn1=(int *) calloc (irnmems,sizeof(int));
//     int *jcn1=(int *) calloc (irnmems,sizeof(int));

    lj2=0;
    for(i=0; i<bivirowsize; i++) {
      bivinzrow1[i]/=bbrowij;
      bivinzrow1[i]++;
    }
    for(i=0; i<bivicolsize; i++) {
      bivinzcol1[i]++;
    }
    //printf("HHHEERE nz %d j %d vs %ld\n",nz,j,vecbivisize);
    for(li=0; li<vecbivisize; ++li) {
      if(vecbivi[li]!=0) {
        vecbivi[lj2]=vecbivi[li];
        irn1[lj2]=bivinzrow1[li/bivicolsize];
        jcn1[lj2]=bivinzcol1[li%bivicolsize];
        lj2++;
      }
    }
    //printf("HHHEERE\n");
//     free(bivinzrow1);
//     bivinzrow1=NULL;
//     //bivinzrow1=realloc (bivinzrow1,1*sizeof(long int));
//     free(bivinzcol1);
//     bivinzcol1=NULL;
    //bivinzcol1=realloc (bivinzcol1,1*sizeof(long int));
    nz0=lj2;
    nrow=submatAij[j4]->rmap->n;
    if(nrow>nrowmax)nrowmax=nrow;
    ncol=submatAij[j4]->cmap->n;
    if(ncol>ncolmax)ncolmax=ncol;
    printf("D ncol %d nrow %d nz %d\n",ncol,nrow,nz);
    for(i=1; i<nrow; i++) {
      for(j=ai[i-1]; j<ai[i]; j++) {
        vecbivi[lj2]=vals[j];
        irn1[lj2]=i;
        jcn1[lj2]=aj[j]+1;
        lj2++;
      }
    }
    for(j=ai[i-1]; j<nz; j++) {
      vecbivi[lj2]=vals[j];
      irn1[lj2]=i;
      jcn1[lj2]=aj[j]+1;
      lj2++;
    }
    //printf("in C %ld\n",biviindx0[1]);
    //vecbivi=realloc(vecbivi,lj2*sizeof(ha_cgetype));
    //if(vecbivisize<ldsize)vecbivisize=ldsize;
    nz1=lj2;
    if(nz1>nz1max)nz1max=nz1;
    MatDestroy(&submatAij[j4]);//1
    //begin save for sol
    if(j3<10)strcpy(j1name,"000");
    if(j3<100&&j3>9)strcpy(j1name,"00");
    if(j3<1000&&j3>99)strcpy(j1name,"0");
    if(j3>=1000)j1name[0]='\0';
    sprintf(filename, "%d",j3);
    strcat(j1name,filename);
    if(isLinux==0) {
      strcpy(filename,"_bivi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fwrite(vecbivi, sizeof(ha_cgetype),nz1, presolfile);
      if(fwrt==0)printf("File write error.");
      fclose(presolfile);
      strcpy(filename,"_rbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fwrite(irn1, sizeof(int), nz1, presolfile);
      if(fwrt== 0) printf("File write error.");
      fclose(presolfile);
      strcpy(filename,"_cbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fwrite(jcn1, sizeof(int), nz1, presolfile);
      //if(rank==2)for(i=0;i<ncol;i++)printf("rank %d j1 %d col %d\n",rank,j3,jcn1[i]);
      if(fwrt== 0) printf("File write error.");
      fclose(presolfile);
    }
    else {
      strcpy(filename,"_bivi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      fd1 = open(filename,O_WRONLY|O_CREAT|O_TRUNC, 0644);
      fwrt = write(fd1,vecbivi, sizeof(ha_cgetype)*nz1);
      close(fd1);
      strcpy(filename,"_rbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      fd1 = open(filename,O_WRONLY|O_CREAT|O_TRUNC, 0644);
      fwrt = write(fd1,irn1, sizeof(int)*nz1);
      close(fd1);
      strcpy(filename,"_cbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      fd1 = open(filename,O_WRONLY|O_CREAT|O_TRUNC, 0644);
      fwrt = write(fd1,jcn1, sizeof(int)*nz1);
      close(fd1);
    }
    free(irn1);//1
    irn1=NULL;
    free(jcn1);//1
    jcn1=NULL;
    }
    //end save for sol

    nz=nz1;
    //obiviindx1=realloc(obiviindx1,sizeof(long int));
    //obiviindx0=realloc(obiviindx0,sizeof(long int));
//     ldsize=ceil((laDi/100.0)*nz1);//second estimate
//     if(irnmems>ldsize){//just in case first and second estimate different
//       irn1=realloc(irn1,ldsize*sizeof(int));
//       jcn1=realloc(jcn1,ldsize*sizeof(int));
//     }else{
//       if(irnmems<ldsize){
//       if(irnmems>nz1){
//         irn1=realloc(irn1,nz1*sizeof(int));
//         jcn1=realloc(jcn1,nz1*sizeof(int));
//       }
//       irn1=realloc(irn1,ldsize*sizeof(int));
//       jcn1=realloc(jcn1,ldsize*sizeof(int));
//       }
//     }
//     if(vecbivisize>ldsize)vecbivi=realloc(vecbivi,ldsize*sizeof(ha_cgetype));
//     else{
//       if(vecbivisize<ldsize){
//       if(lj2<vecbivisize)vecbivi=realloc(vecbivi,lj2*sizeof(ha_cgetype));
//       vecbivi=realloc(vecbivi,ldsize*sizeof(ha_cgetype));
//       }
//     }
    insize[j4*insizes]=nrow;
    insize[j4*insizes+1]=ncol;
    insize[j4*insizes+2]=nz;
    insize[j4*insizes+10]=rank;
    insize[j4*insizes+11]=j4;
    insize[j4*insizes+9]=laDi;
    insize[j4*insizes+15]=7+windx;//7+rank+jthrd+windx;//(j4+jthrd)%90+7;
//     windx++;
//     if(windx==eindx)windx=bindx;
    //for(i=0;i<nz;i++)printf("aii %f irn1 %d jcn %d\n",vecbivi[i],irn1[i],jcn[i]);
    //printf("ok!!!!!!!! vecbiuisize %d nrow %d\n",ncol,nrow);
      //int *insized=(int *) calloc (5+nreg*insizes,sizeof(int));
      insized[0]=nrow;
      insized[1]=ncol;
      insized[2]=nz;
      insized[4]=laDi;
      //if(rank==2)for(i=0;i<ncol;i++)printf("arank %d j1 %d col %d sct %d svb %d\n",rank,j3,jcn[i],(int)sizeof(PetscReal),(int)sizeof(ha_cgetype));
//       if(cntl6==0&&SORD==0)cntl6in=0.3;
//       else cntl6in=cntl6;
//       spec51_rank_(insized,&cntl6in,irn1,jcn1,vecbivi);
      //insized[3]=(int)(0.9*insized[3]);
      //insized=realloc(insized,(5+nreg*insizes)*sizeof(int));
      for(i=0; i<nreg; i++) {
        for(j=0; j<insizes; j++) {
          insized[i*insizes+5+j]=insize[(i+j3*(nreg+1))*insizes+j];
        }
      }
      if(j3<10)strcpy(j1name,"000");
      if(j3<100&&j3>9)strcpy(j1name,"00");
      if(j3<1000&&j3>99)strcpy(j1name,"0");
      if(j3>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j3);
      strcat(j1name,filename);
      strcpy(filename,"_rank");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      printf("Rank %d RRRRR! ncol %d nrow %d rank %d nz %d cntl6 %lf\n",rank,nrow,ncol,insized[3],nz,cntl6);
      fwrt=fwrite(insized, sizeof(int), 5+nreg*insizes, presolfile);
      if(fwrt==0)printf("File write error.");
      fclose(presolfile);
//       strcpy(filename,"_row");
//       strcat(filename,rankname);
//       strcat(filename,j1name);
//       strcat(filename,".bin");
//       if((presolfile=fopen(filename, "wb"))==NULL) {
//         printf("Cannot open file.\n");
//       }
//       fwrt=fwrite(irn1, sizeof(int), nrow, presolfile);
//       if(fwrt== 0) printf("File write error.");
//       fclose(presolfile);
//       strcpy(filename,"_col");
//       strcat(filename,rankname);
//       strcat(filename,j1name);
//       strcat(filename,".bin");
//       if((presolfile=fopen(filename, "wb"))==NULL) {
//         printf("Cannot open file.\n");
//       }
//       fwrt=fwrite(jcn1, sizeof(int), ncol, presolfile);
//       //if(rank==2)for(i=0;i<ncol;i++)printf("rank %d j1 %d col %d\n",rank,j3,jcn1[i]);
//       if(fwrt== 0) printf("File write error.");
//       fclose(presolfile);
      //free(insized);//1
    insize[j4*insizes+15]=0;
    //printf("ok!!!!!!!!\n");
    //vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
  }
      free(irn);//1
      irn=NULL;
      free(jcn);//1
      jcn=NULL;
      free(values);//1
    free(jcnb1);
    jcnb1=NULL;
    free(sol48);
    sol48=NULL;
    free(b48);
    b48=NULL;
    free(keep);
    keep=NULL;
    free(w51);
    w51=NULL;
    free(iw51);
    iw51=NULL;
    free(insized);
    free(ccolnorms);//1
    ccolnorms=NULL;
    free(bivinzcol);//1
    bivinzcol=NULL;
    free(bivinzrow);//1
    bivinzrow=NULL;
    free(bivinzrow1);//1
    bivinzrow1=NULL;
    free(bivinzcol1);//1
    bivinzcol1=NULL;
    free(vecbivi);//1
    vecbivi=NULL;
  //free(nfwrt);
  ldsize=ceil((laDi/100.0)*nz1max);
  ldsize+=10;
  irn1=realloc(irn1,ldsize*sizeof(int));
  if(irn1==NULL)printf("Error!!!! No memory left for irn1!\n");
  jcn1=realloc(jcn1,ldsize*sizeof(int));
  if(jcn1==NULL)printf("Error!!!! No memory left for jcn1!\n");
  irn=realloc(irn,nrowmax*sizeof(int));
  if(irn==NULL)printf("Error!!!! No memory left for irn!\n");
  jcn=realloc(jcn,ncolmax*sizeof(int));
  if(jcn==NULL)printf("Error!!!! No memory left for jcn!\n");
  keep=(int *) malloc ((nrowmax+5*ncolmax+4*ncolmax+7)*sizeof(int));//ICNTL(6)=1
  if(keep==NULL)printf("Error!!!! No memory left for keep!\n");
  //ha_cgetype *w51;
  if(nrowmax>ncolmax) w51=(ha_cgetype *) malloc ((5*nrowmax)*sizeof(ha_cgetype));
  else w51=(ha_cgetype *) malloc ((5*ncolmax)*sizeof(ha_cgetype));
  if(w51==NULL)printf("Error!!!! No memory left for w51!\n");
  iw51=(int *) malloc ((6*nrowmax+3*ncolmax)*sizeof(int));
  if(iw51==NULL)printf("Error!!!! No memory left for iw51!\n");
  vecbivi=realloc(vecbivi,ldsize*sizeof(ha_cgetype));
  if(vecbivi==NULL)printf("Error!!!! No memory left for vecbivi!\n");
  int *insizeda=(int *) calloc (5+nreg*insizes,sizeof(int));
  #pragma omp for schedule (static)
  for(j3=0; j3<nmatint; j3++) {
    j4=nreg+j3*(nreg+1);
    if(j3<10)strcpy(j1name,"000");
    if(j3<100&&j3>9)strcpy(j1name,"00");
    if(j3<1000&&j3>99)strcpy(j1name,"0");
    if(j3>=1000)j1name[0]='\0';
    sprintf(filename, "%d",j3);
    strcat(j1name,filename);
    strcpy(filename,"_rank");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    //printf("filename %s\n",filename);
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Cannot open file. Please run presol by setting -presol 1!\n");
    }
    fwrt=fread(insizeda, sizeof(int), 5+nreg*insizes, presolfile);
    fclose(presolfile);
    remove(filename);
    if( fwrt== 0)printf("File read error. No presol!\n");
    nrow=insizeda[0];
    ncol=insizeda[1];
    nz1=insizeda[2];
    nz=nz1;
      strcpy(filename,"_bivi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fread(vecbivi, sizeof(ha_cgetype), insizeda[2], presolfile);
      if(fwrt== 0) printf("File read error %ld %d.\n",fwrt,insizeda[2]);
      fclose(presolfile);
      strcpy(filename,"_rbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fread(irn1, sizeof(int), insizeda[2], presolfile);
      if(fwrt== 0) printf("File read error.\n");
      fclose(presolfile);
      strcpy(filename,"_cbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fread(jcn1, sizeof(int), insizeda[2], presolfile);
      if(fwrt== 0) printf("File read error.");
      fclose(presolfile);

    insize[j4*insizes]=nrow;
    insize[j4*insizes+1]=ncol;
    insize[j4*insizes+2]=nz;
    insize[j4*insizes+10]=rank;
    insize[j4*insizes+11]=j4;
    insize[j4*insizes+9]=laDi;
    insize[j4*insizes+15]=7+windx;//7+rank+jthrd+windx;//(j4+jthrd)%90+7;
    windx++;
    if(windx==eindx)windx=bindx;
      int *insized=(int *) calloc (6,sizeof(int));
      insized[0]=nrow;
      insized[1]=ncol;
      insized[2]=nz;
      insized[4]=laDi;
    ldsize=ceil((laDi/100.0)*nz);
      insized[5]=ldsize;
      //if(rank==2)for(i=0;i<ncol;i++)printf("arank %d j1 %d col %d sct %d svb %d\n",rank,j3,jcn[i],(int)sizeof(PetscReal),(int)sizeof(ha_cgetype));
      if(cntl6==0&&SORD==0)cntl6in=0.3;
      else cntl6in=cntl6;
      spec51m_rank_(insized,&cntl6in,irn1,jcn1,vecbivi,irn,jcn,keep,w51,iw51);
      //insized[3]=(int)(0.9*insized[3]);
      insized=realloc(insized,(5+nreg*insizes)*sizeof(int));
      for(i=0; i<nreg; i++) {
        for(j=0; j<insizes; j++) {
          insized[i*insizes+5+j]=insize[(i+j3*(nreg+1))*insizes+j];
        }
      }
      if(j3<10)strcpy(j1name,"000");
      if(j3<100&&j3>9)strcpy(j1name,"00");
      if(j3<1000&&j3>99)strcpy(j1name,"0");
      if(j3>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j3);
      strcat(j1name,filename);
      strcpy(filename,"_rank");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      printf("Rank %d RRRRR! ncol %d nrow %d rank %d nz %d cntl6in %lf\n",rank,nrow,ncol,insized[3],nz,cntl6in);
      fwrt=fwrite(insized, sizeof(int), 5+nreg*insizes, presolfile);
      if(fwrt==0)printf("File write error.");
      fclose(presolfile);
      strcpy(filename,"_row");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fwrite(irn, sizeof(int), nrow, presolfile);
      if(fwrt== 0) printf("File write error.");
      fclose(presolfile);
      strcpy(filename,"_col");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fwrite(jcn, sizeof(int), ncol, presolfile);
      //if(rank==2)for(i=0;i<ncol;i++)printf("rank %d j1 %d col %d\n",rank,j3,jcn1[i]);
      if(fwrt== 0) printf("File write error.");
      fclose(presolfile);
      free(insized);//1
    insize[j4*insizes+15]=0;
    //printf("ok!!!!!!!!\n");
    //vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
//     free(vecbivi);//1
//     free(irn1);//1
//     irn1=NULL;
//     free(jcn1);//1
//     jcn1=NULL;
//     free(ccolnorms);//1
//     ccolnorms=NULL;
//     free(bivinzcol);//1
//     bivinzcol=NULL;
//     free(bivinzrow);//1
//     bivinzrow=NULL;
//     free(bivinzrow1);//1
//     bivinzrow1=NULL;
//     free(bivinzcol1);//1
//     bivinzcol1=NULL;
  }
    free(vecbivi);//1
    free(irn1);//1
    free(jcn1);//1
    free(irn);//1
    free(jcn);//1
    free(keep);
    free(w51);
    free(iw51);
  free(insizeda);
  }
  //omp_set_num_threads(mymaxnumthrd);
  //if(presol)return 0;
  //free(npresolfilea);
  //free(npresolfilea1);
  //free(npresolfileb);
  //free(npresolfileb1);
  //free(npresolfilec);
  //free(npresolfilec1);
  ierr = PetscFree(submatAij);//1
  printf("Rank %d End of preparation!\n",rank);
  MPI_Barrier(PETSC_COMM_WORLD);

  //bbrowij=submatB[0]->rmap->n;
  //printf("j3 %d bbrowij %d\n",j3,bbrowij);
  //ha_cgetype *vecbivi= (ha_cgetype *) calloc (vecbivisize,sizeof(ha_cgetype));
//  for(j1=0; j1<nmatint; j1++) {
//    MatDestroy(&submatCT);
//    free(yi1[j1]);//1
//    MatDestroy(&submatB1[j1][0]);//submatBT);//1
//    ierr = PetscFree(submatB1[j1]);//1
//    MatDestroy(&submatB2[j1][0]);//submatBT);//1
//    ierr = PetscFree(submatB2[j1]);//1
  //ierr = PetscFree(ai12);
  //ierr = PetscFree(aj12);
  //ierr = PetscFree(vals12);
//  }
//  ierr = PetscFree(submatB1);//1
//  ierr = PetscFree(submatB2);//1
//  CHKERRQ(ierr);
//  free(yi1);//1
  //vecbivi=realloc(vecbivi,j2*sizeof(ha_cgetype));
//  if(rank==mpisize-1) {
//    MatDestroy(&submatD[0]);//1
//    //free(irn1);
//    //free(jcn);
//    //free(insizeD);
//    //x0=(ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));//realloc (x0,*sizeof(ha_cgetype));
//  }
//  free(vecbiui);//1
  //free(obiviindx1);//1
  //free(obiviindx0);//1
  //free(vecbivi);//1
//  free(indicesC);//1
//  PetscFree(submatD);//1
//  for(j1=0; j1<nmatint; j1++) {
//    for(i=0; i<nreg; i++) {
//      MatDestroy(&submatCij[j1*(nreg+1)+i][0]);//1
//      MatDestroy(&submatBij[j1*(nreg+1)+i][0]);//1
//    }
//  }
  for(i=0; i<nmatin; i++) {
    ierr = PetscFree(submatCij[i]);//1
    ierr = PetscFree(submatBij[i]);//1
  }
//  ierr = PetscFree(colindicesbc1);//1
//  ierr = PetscFree(colindicesbc2);//1
  free(submatCij);//1
  free(submatBij);//1

//  free(xi1);//1
  free(insize);//1
  //ierr = PetscFree(submatC);//1
  //free(begblock);//1
  //timestr=clock();
  //free(cirank);//1
  //free(birank);//1
  //free(aord);//1
  //free(ndblockinrank);//1
  return 0;
}

int HaNDBBDParPre1(Mat A, Vec b, ha_cgetype *x1, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, uvadd ndblock,uvadd nreg,uvadd ntime,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,uvdim laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol) {//,bool iter
  IS *rowindices,*colindices,*rowindicesbc,*colindicesbc1,*colindicesbc2,*Cindices,*Bindices,*BBindices,*BBiindices,*rowBBij,*colBBij;//,*colindicesbcpm,*colindicesbcpm1
  const PetscInt *nindices;
  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nmatint,nmatinplust,nmatminust,nrowcolin,i,i1,j,j0,j1,j2,j3,j4,j5,j6,j7,l0,l1,l2,l3,l4,l5,rank1,proc1=0,nnzmax,j1nz,j1irnbs;//,sumrowcolin
  Mat *submatAij,*submatC,*submatD,submatBT,submatCT,*submatBB,*submatBBij;//,*submatCij,*submatBij;,*submatB
  Vec vecxd;
  PetscInt *ai,*aj,*a1i,*a1j,*a2i,*a2j,*aic,*ajc;
  PetscReal cntl6in;
  uvadd lasize,ldsize;
  long int *biviindx1,*biviindx0,j7l;
  long int *obiviindx1= (long int*)calloc(1,sizeof(long int));
  long int *obiviindx0= (long int*)calloc(1,sizeof(long int));
  size_t freadresult,fwrt;
  obiviindx1[0]=-1;
  obiviindx0[0]=-1;
  PetscInt nrow,ncol,nz,nrowc,ncolc,nrowb,ncolb,nzc,nzv;
  MPI_Status   status;
  clock_t timestr,timeend;
  int j2int,la1;
  bool ifremove=false;
  PetscScalar *vals,*val1s,*val2s,*valsc,vecval;
  ha_cgetype vval;
  PetscErrorCode ierr;
  PetscViewer viewer;
  FILE* fp1,*fp2,*fp3,*presolfile;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  int *ndblockinrank= (int *) calloc (mpisize,sizeof(int));
  char filename[1024],rankname[1024],j1name[1024],tempchar[1024];//,fn1[1024],fn2[1024],fn3[1024]
  if(rank<10)strcpy(rankname,"000");
  if(rank<100&&rank>9)strcpy(rankname,"00");
  if(rank<1000&&rank>99)strcpy(rankname,"0");
  if(rank>=1000)rankname[0]='\0';
  sprintf(filename, "%d",rank);
  strcat(rankname,filename);
  int nmatinBB=1;
  nmatint=(uvadd)ntime/mpisize;
  nmatminust=nmatint;
  if(nmatint*mpisize<ntime)nmatinplust=nmatint+1;
  else nmatinplust=nmatint;
  for(i=0; i<mpisize; i++)if(rank+1<=ntime-mpisize*nmatint)nmatint++;
  nmatin=(nreg+1)*nmatint;
  nmatinplus=(nreg+1)*nmatinplust;
  begblock[rank]=nmatin;
  printf("rank %d nmatin %d nmatint %d nmplus %d\n",rank,nmatin,nmatint,nmatinplus);
  for(i=0; i<mpisize; i++) {
    j=rank;
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    j=j2;
    MPI_Bcast(&begblock[j2],1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  j2=begblock[0];
  begblock[0]=0;
  for(i=1; i<mpisize; i++) {
    j1=begblock[i];
    begblock[i]=begblock[i-1]+j2;
    j2=j1;
  }
  //for(i=0; i<mpisize; i++)printf("rank %d beg block %d\n",rank,begblock[i]);
  j2=nmatin;
  MPI_Bcast(&j2,1, MPI_INT,0, PETSC_COMM_WORLD);
  ndblockinrank[0]=j2;
  for(i=1; i<mpisize; i++) {
    j2=nmatin;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    ndblockinrank[i]=ndblockinrank[i-1]+j2;
  }
  int *birank= (int *) calloc (nmatin,sizeof(int));
  int *cirank= (int *) calloc (nmatin,sizeof(int));
  //for(i=0; i<nmatin; i++)cirank[i]=-1;
  int *aord= (int *) calloc (nmatin,sizeof(int));
  for(i=0; i<ndblock; i++) {
    j2=i-ndblockinrank[rank];
    if(j2<0&&j2>=-nmatin) {
      j1=i/(nreg+1);
      j=i-(nreg+1)*j1;
      if(rank==0)j3=0;
      else j3=ndblockinrank[rank-1];
      if(j<nreg)cirank[j2+nmatin]=(nreg+1)*j1+nreg-j3;
      else cirank[j2+nmatin]=-1;
    }
  }
  j1=0;
  for(j=0; j<mpisize; j++)if(j<rank)j1+=ndblockinrank[j];
  for(i=0; i<nmatin; i++) {
    j2=(j1+i)/(nreg+1);
    if(i-j2*(nreg+1)!=0)aord[i]=1;
  }
  ierr = PetscMalloc(nmatin*sizeof(IS **),&rowindices);
  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*sizeof(IS **),&rowindicesbc);
//  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&colindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbc1);
  CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbc2);
  //CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*sizeof(IS **),&Cindices);
//  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*sizeof(IS **),&Bindices);
//  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatinBB*sizeof(IS **),&BBindices);
//  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatAij);
  CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatint*sizeof(Mat *),&submatC);
  //CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatinBB*sizeof(Mat *),&submatBB);
//  CHKERRQ(ierr);
  Mat** submatBij = (Mat**)calloc(nmatin,sizeof(Mat*));
  for (i=0; i<nmatin; i++) {
    submatBij[i] = (Mat*)calloc(1,sizeof(Mat));
  }
  Mat** submatCij = (Mat**)calloc(nmatin,sizeof(Mat*));
  for (i=0; i<nmatin; i++) {
    submatCij[i] = (Mat*)calloc(1,sizeof(Mat));
  }
//  Mat** submatB1 = (Mat**)calloc(nmatint,sizeof(Mat*));
//  for (i=0; i<nmatint; i++) {
//    submatB1[i] = (Mat*)calloc(1,sizeof(Mat));
//  }
//  Mat** submatB2 = (Mat**)calloc(nmatint,sizeof(Mat*));
//  for (i=0; i<nmatint; i++) {
//    submatB2[i] = (Mat*)calloc(1,sizeof(Mat));
//  }

//  ha_cgetype **yi2= (ha_cgetype**)calloc(nmatin,sizeof(ha_cgetype*));
//  for (i=0; i<nmatin; i++) {
//    yi2[i] = (ha_cgetype*)calloc(ha_ndblocks[i+begblock[rank]],sizeof(ha_cgetype));
//  }
//  sumrowcolin=0;
//  for(j=0; j<ndblock; j++)sumrowcolin+=ha_ndblocks[j];
//  PetscInt *offblock= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
//  PetscInt *offblockrow= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
//  offblock[0]=0;
//  for(j=0; j<ndblock; j++)offblock[j+1]=offblock[j]+countvarintra1[j+1]-countvarintra1[j]-ha_ndblocks[j];
//  offblockrow[0]=0;
//  for(j=0; j<ndblock; j++)offblockrow[j+1]=offblockrow[j]+counteqnoadd[j]-ha_ndblocks[j];

  char **fn01= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn01[i] = (char*)calloc(1024,sizeof(char));
  char **fn02= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn02[i] = (char*)calloc(1024,sizeof(char));
  char **fn03= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn03[i] = (char*)calloc(1024,sizeof(char));

  timestr=clock();
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      if(j1<10)strcpy(j1name,"000");
      if(j1<100&&j1>9)strcpy(j1name,"00");
      if(j1<1000&&j1>99)strcpy(j1name,"0");
      if(j1>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j1);
      strcat(j1name,filename);
      strcpy(filename,"_irnv");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn01[j1],filename);
      strcpy(filename,"_keep");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn02[j1],filename);
      strcpy(filename,"_vav");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn03[j1],filename);
    }
  }
  //printf("OK012345 rank %d mpis %d sizeof int %ld size of long int %ld\n",rank,mpisize,sizeof(int),sizeof(long int));
//  for(j=0; j<mpisize; j++) {
//    //printf("j %d\n",j);
//    for(j3=0; j3<nmatinplus; j3++) {
//      if(j3>=nmatin)i=nmatin-1;
//      else i=j3;
//      //printf("rank %d OK1!!! i %d j %d nmatin %d\n",rank,i,j,nmatin);
//      bfirst=counteq[i+begblock[j]];
//      bend=ha_ndblocks[i+begblock[j]];
//      MPI_Bcast(&bfirst,1, MPI_INT,j, PETSC_COMM_WORLD);
//      MPI_Bcast(&bend,1, MPI_INT,j, PETSC_COMM_WORLD);
//      PetscInt *indices= (PetscInt *) calloc (bend,sizeof(PetscInt));
//      ha_cgetype *yi0 = (ha_cgetype*)calloc(bend,sizeof(ha_cgetype));
//      for(j1=0; j1<bend; j1++) { //ha_ndblocks[i+begblock[j]]
//        j2=ha_rows[bfirst+j1];//+bfirst;
//        if(j2>=Istart&&j2<Iend)indices[j1]=j2;
//        else indices[j1]=-1;
//      }
//      //printf("rank %d i %d j %d bf %d be %d\n",rank,i,j,bfirst,bend);
//      VecGetValues(b,bend,indices,yi0);
//      HaReduce(yi0,(forint)bend,mpisize,rank,j);
//      if(rank==j&&i<nmatin) {
//        memcpy (yi2[i],yi0,bend*sizeof(ha_cgetype));
//      }
//      free(indices);//1
//      free(yi0);//1
//    }
//  }
//  ha_cgetype **yi1= (ha_cgetype**)calloc(nmatint,sizeof(ha_cgetype*));
//  for (i=0; i<nmatint; i++) {
//    j1=0;
//    for (j=0; j<nreg+1; j++) {
//      j1+=ha_ndblocks[i*(nreg+1)+j+begblock[rank]];
//    }
//    yi1[i] = (ha_cgetype*)calloc(j1,sizeof(ha_cgetype));
//  }
//  for (i=0; i<nmatint; i++) {
//    j1=0;
//    for (j=0; j<nreg+1; j++) {
//      memcpy (&yi1[i][j1],yi2[i*(nreg+1)+j],ha_ndblocks[i*(nreg+1)+j+begblock[rank]]*sizeof(ha_cgetype));
//      free(yi2[i*(nreg+1)+j]);//1
//      j1+=ha_ndblocks[i*(nreg+1)+j+begblock[rank]];
//    }
//  }
//  free(yi2);//1
  //printf("OK12345 rank %d\n",rank);
  j1=0;
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (ha_ndblocks[i+begblock[rank]],sizeof(PetscInt));
    //printf("ha_row 0 %d 588 %d\n",ha_rows[0],ha_rows[588]);
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_rows[bfirst+j];//+bfirst;
    }

    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
    bfirst=countvarintra1[i+begblock[rank]];
    //if(rank==0)printf("rank %d i %d bfirst %d block %d\n",rank,i,bfirst,ha_ndblocks[i+begblock[rank]]);
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_cols[bfirst+j];//+bfirst;
    }
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,colindices+i);
    free(indices);//1
  }
  for(i=0; i<nmatint; i++) {
    j2=0;
    for(j1=0; j1<nreg+1; j1++)j2+=ha_ndblocks[i*(nreg+1)+j1+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (j2,sizeof(PetscInt));
//    j2=0;
//    for(j1=0; j1<nreg+1; j1++) {
//      bfirst=counteq[i*(nreg+1)+j1+begblock[rank]];
//      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
//        indices[j2]=ha_rows[bfirst+j];//+bfirst;
//        j2++;
//      }
//    }
//    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,rowindicesbc+i);
    j2=0;
    for(j1=0; j1<nreg; j1++) {
      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
        indices[j2]=ha_cols[bfirst+j];//+bfirst;
        j2++;
      }
    }
    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbc1+i);
//    j2=0;
//    for(j1=nreg; j1<nreg+1; j1++) {
//      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
//      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
//        indices[j2]=ha_cols[bfirst+j];//+bfirst;
//        j2++;
//      }
//    }
    //ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbc2+i);
    free(indices);
  }
//  printf("0okkkk\n");
//  PetscInt *indicesC= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
//  PetscInt *indicesB= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
//  for(i=0; i<ndblock; i++) {
//    bend=countvarintra1[i+1];
//    for(j=offblock[i]; j<offblock[i+1]; j++) {
//      j1=j-offblock[i];
//      indicesC[j]=ha_cols[bend-j1-1];//+bfirst;
//    }
//  }
//  j1=offblock[ndblock];
//  for(i=countvarintra1[ndblock]; i<VecSize; i++) {
//    indicesC[j1]=i;
//    j1++;
//  }
//  printf("1okkkk\n");
//  for(i=0; i<ndblock; i++) {
//    bend=counteq[i]+counteqnoadd[i];
//    for(j=offblockrow[i]; j<offblockrow[i+1]; j++) {
//      j1=j-offblockrow[i];
//      indicesB[j]=ha_rows[bend-j1-1];//+bfirst;
//    }
//  }
//  j1=offblockrow[ndblock];
//  for(i=counteq[ndblock-1]+counteqnoadd[ndblock-1]; i<VecSize; i++) {
//    indicesB[j1]=i;
//    j1++;
//  }
//  printf("2okkkk\n");
//  printf("rank1 %d\n",rank);
//  for(i=0; i<VecSize-sumrowcolin; i++)if(indicesB[i]>=Istart) {
//      j=i;
//      break;
//    }
//  for(j2=i; j2<VecSize-sumrowcolin; j2++)if(indicesB[j2]>=Iend) {
//      j1=j2;
//      break;
//    }
//  if(j2==VecSize-sumrowcolin)j1=VecSize-sumrowcolin;
//  printf("rank %d j %d j1 %d istart %d iend %d\n",rank,j,j1,Istart,Iend);
//  ha_cgetype *vecbiui= (ha_cgetype *) calloc (VecSize-sumrowcolin,sizeof(ha_cgetype));
//  VecGetValues(b,j1-j,&indicesB[j],&vecbiui[j]);//implicite yd save mem
  //ierr = VecDestroy(&b);
  //CHKERRQ(ierr);
//  printf("rank2 %d\n",rank);
//  for(i=0; i<nmatint; i++) {
//    ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesC,PETSC_COPY_VALUES,Cindices+i);
//    ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesB,PETSC_COPY_VALUES,Bindices+i);
//  }
//  free(indicesB);//1
//  free(offblock);//1
//  free(offblockrow);//1
//  PetscInt *indexBB= (PetscInt *) calloc (VecSize,sizeof(PetscInt));
//  for(i=0; i<VecSize; i++)indexBB[i]=i;
//  ISCreateGeneral(PETSC_COMM_SELF,VecSize,indexBB,PETSC_COPY_VALUES,BBindices);
  printf("Begin preparation rank %d\n",rank);
  ierr = MatCreateSubMatrices(A,nmatin,rowindices,colindices,MAT_INITIAL_MATRIX,&submatAij);
  //if(rank==0)printf("rank %d i %d bbrowij %d\n",rank,(nmatint-1)*(nreg+1)+nreg,submatAij[(nmatint-1)*(nreg+1)+nreg]->cmap->n);
  int insizes=17;
  int *insize=(int *) calloc (insizes*nmatin,sizeof(int));
  printf("rank %d A\n",rank);
  PetscInt *indicesbbij=(int *) calloc (1,sizeof(PetscInt));
  ierr = PetscMalloc(1*sizeof(IS **),&rowBBij);
  CHKERRQ(ierr);
  ierr = PetscMalloc(1*sizeof(IS **),&colBBij);
  CHKERRQ(ierr);

  for(i=0; i<nmatminust; i++) {
    j2=i*(nreg+1)+nreg;
    ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindicesbc1+i,MAT_INITIAL_MATRIX,&submatBBij);
    nrow=submatBBij[0]->rmap->n;
    ncol=submatBBij[0]->cmap->n;
    if(nrow>ncol) {
      indicesbbij=realloc(indicesbbij,nrow*sizeof(PetscInt));
      for(j=0; j<nrow; j++)indicesbbij[j]=j;
    }
    else {
      indicesbbij=realloc(indicesbbij,ncol*sizeof(PetscInt));
      for(j=0; j<ncol; j++)indicesbbij[j]=j;
    }
    ISCreateGeneral(PETSC_COMM_SELF,nrow,indicesbbij,PETSC_COPY_VALUES,rowBBij);
    j3=0;
    for(j=0; j<nreg; j++) {
      j1=i*(nreg+1)+j;
      ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
      ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[j1+begblock[rank]],indicesbbij+j3,PETSC_COPY_VALUES,colBBij);
      ierr = MatCreateSubMatrices(submatBBij[0],1,rowBBij,colBBij,MAT_INITIAL_MATRIX,&submatBij[j1]);
      ierr = ISDestroy(&colBBij[0]);
      j3+=ha_ndblocks[j1+begblock[rank]];
      //printf("rank %d i %d j3 %d bbrowij %d\n",rank,i,j,submatCij[j1][0]->cmap->n);
    }
    ierr = MatDestroy(&submatBBij[0]);
    ierr = ISDestroy(&rowBBij[0]);
  }

//  for(i=0; i<nmatminust; i++) {
//      j2=i*(nreg+1)+nreg;
//    for(j=0; j<nreg; j++) {
//      j1=i*(nreg+1)+j;
////      ierr = ISSort(rowindices[j2]);
////      ierr = ISSort(rowindices[j1]);
////      ierr = ISSort(colindices[j1]);
////      ierr = ISSort(colindices[j2]);
//      ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
//      //printf("i %d j %d C\n",i,j);
//      ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindices+j1,MAT_INITIAL_MATRIX,&submatBij[j1]);
//      //printf("i %d j %d B\n",i,j);
//      //printf("rank %d i %d j3 %d bbrowij %d\n",rank,i,j,submatCij[j1][0]->cmap->n);
//    }
//  }
//  //printf("rank %d A1\n",rank);
  if(nmatint<nmatinplust) for(i=0; i<nreg; i++) {//to create again! MatCreateSubMatrices must be the same across
      MatDestroy(&submatCij[(nmatint-1)*(nreg+1)+i][0]);
      MatDestroy(&submatBij[(nmatint-1)*(nreg+1)+i][0]);
    }
  j2=(nmatint-1)*(nreg+1)+nreg;
  ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindicesbc1+(nmatint-1),MAT_INITIAL_MATRIX,&submatBBij);
  nrow=submatBBij[0]->rmap->n;
  ncol=submatBBij[0]->cmap->n;
  if(nrow>ncol) {
    indicesbbij=realloc(indicesbbij,nrow*sizeof(PetscInt));
    for(j=0; j<nrow; j++)indicesbbij[j]=j;
  }
  else {
    indicesbbij=realloc(indicesbbij,ncol*sizeof(PetscInt));
    for(j=0; j<ncol; j++)indicesbbij[j]=j;
  }
  ISCreateGeneral(PETSC_COMM_SELF,nrow,indicesbbij,PETSC_COPY_VALUES,rowBBij);
  j3=0;
  for(j=0; j<nreg; j++) {
    j1=(nmatint-1)*(nreg+1)+j;
    ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[j1+begblock[rank]],indicesbbij+j3,PETSC_COPY_VALUES,colBBij);
    ierr = MatCreateSubMatrices(submatBBij[0],1,rowBBij,colBBij,MAT_INITIAL_MATRIX,&submatBij[j1]);
    ierr = ISDestroy(&colBBij[0]);
    j3+=ha_ndblocks[j1+begblock[rank]];
    //printf("rank %d i %d j3 %d bbrowij %d\n",rank,i,j,submatCij[j1][0]->cmap->n);
  }
  ierr = MatDestroy(&submatBBij[0]);
  ierr = ISDestroy(&rowBBij[0]);
  free(indicesbbij);
  PetscFree(rowBBij);
  PetscFree(colBBij);

//  for(j=0; j<nreg; j++) {
//    j1=(nmatint-1)*(nreg+1)+j;
//    j2=(nmatint-1)*(nreg+1)+nreg;
//    ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
//    ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindices+j1,MAT_INITIAL_MATRIX,&submatBij[j1]);
//    //printf("rank %d i %d j3 %d bbrowij %d\n",rank,i,j1,submatCij[j1][0]->cmap->n);
//  }
  //printf("rank %d A\n",rank);
  //ierr = MatCreateSubMatrices(A,nmatint,rowindicesbc,Cindices,MAT_INITIAL_MATRIX,&submatC);
  //printf("C\n");
//  ierr = MatCreateSubMatrices(A,nmatinBB,Bindices,BBindices,MAT_INITIAL_MATRIX,&submatBB);
  //ierr = MatDestroy(&A);
  //CHKERRQ(ierr);
//  PetscInt BBrow,BBcol;//,CCrow,CCcol;
//  MatGetSize(submatBB[0],&BBrow,&BBcol);
//  //MatGetSize(submatC[0],&CCrow,&CCcol);
//  printf("rank %d bbrow %d bbcol %d\n",rank,BBrow,BBcol);//,CCrow,CCcol);
//  PetscInt *indexBBi= (PetscInt *) calloc (BBrow,sizeof(PetscInt));
//  for(i=0; i<BBrow; i++)indexBBi[i]=i;
//  for(i=0; i<nmatint; i++) {
//    ierr = ISDestroy(&Bindices[i]);
//    ISCreateGeneral(PETSC_COMM_SELF,BBrow,indexBBi,PETSC_COPY_VALUES,Bindices+i);
//  }
  //ierr = MatCreateSubMatrices(submatBB[0],nmatint,Bindices,colindicesbc,MAT_INITIAL_MATRIX,&submatB);
  //char           ordering[256] = MATORDERINGRCM;
//  for(i=0; i<nmatint; i++) {
//    ierr = MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc1[i],MAT_INITIAL_MATRIX,&submatB1[i][0]);
//    ierr = ISDestroy(&colindicesbc1[i]);
//    ierr = MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc2[i],MAT_INITIAL_MATRIX,&submatB2[i][0]);
//    ierr = ISDestroy(&colindicesbc2[i]);
//  }
//  ierr = PetscFree(colindicesbc1);//1
//  ierr = PetscFree(colindicesbc2);//1
//  printf("B\n");
//  if(rank==mpisize-1)proc1=1;
//  ierr = PetscMalloc(proc1*sizeof(Mat *),&submatD);
//  CHKERRQ(ierr);
//  printf("D\n");
//  MatCreateSubMatrices(submatBB[0],proc1,Bindices,Cindices,MAT_INITIAL_MATRIX,&submatD);
  printf("End partitionning rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&rowindices[i]);//1
    CHKERRQ(ierr);
    ierr = ISDestroy(&colindices[i]);//1
    CHKERRQ(ierr);
  }
//  for (i=0; i<nmatint; i++) {
//    ierr = ISDestroy(&Cindices[i]);//1
//    CHKERRQ(ierr);
//    ierr = ISDestroy(&Bindices[i]);//1
//    CHKERRQ(ierr);
//  }
  PetscFree(colindices);//1
//  ierr = ISDestroy(&BBindices[0]);//1
  PetscFree(rowindices);//1

//  for (i=0; i<nmatint; i++) {
//    ierr = ISDestroy(&rowindicesbc[i]);//1
//    CHKERRQ(ierr);
//  }
//  PetscFree(rowindicesbc);//1

//  PetscFree(Cindices);//1
//  PetscFree(Bindices);//1
//  PetscFree(BBindices);//1
//  PetscFree(indexBB);//1
  //PetscFree(indexBBi);//1
//  ierr = MatDestroy(&submatBB[0]);//1
//  PetscFree(submatBB);//1
//  CHKERRQ(ierr);
//  printf("OK1 rank %d sumrowcolin %d\n",rank,sumrowcolin);
  //ha_cgetype *xi1 = (ha_cgetype*)calloc(sumrowcolin,sizeof(ha_cgetype));
  long int *bivinzrow=(long int *) calloc (1,sizeof(long int));
  PetscInt *bivinzcol=(PetscInt *) calloc (1,sizeof(PetscInt));
  PetscReal *ccolnorms=(PetscReal *) calloc (1,sizeof(PetscReal));
  //ccolnorms=realloc(ccolnorms,1*sizeof(PetscReal));
  int bivirowsize=1,bivicolsize=1,bbrowij,ddrowi;
  long int vecbivisize=0,halfvec;
  ha_cgetype *xi1point;
  long int *bivinzrow1=(long int *) calloc (1,sizeof(long int));
  int *bivinzcol1=(int *) calloc (1,sizeof(int));
  //ha_cgetype *vecbivi= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
  uvadd xi1indx=0;
  int nz0,nz1,nz2,nz3=j,halfj2;

  for(j3=0; j3<nmatint; j3++) {
    bivirowsize=1;
    bivicolsize=1;
    bbrowij=submatBij[j3*(nreg+1)][0]->rmap->n;
    printf("j3 %d bbrowij %d\n",j3,bbrowij);
    ccolnorms=realloc(ccolnorms,bbrowij*sizeof(PetscReal));
    bivinzrow=realloc(bivinzrow,bbrowij*sizeof(long int));
    memset(bivinzrow,0,bbrowij*sizeof(long int));
    bivinzcol=realloc(bivinzcol,bbrowij*sizeof(PetscInt));
    memset(bivinzcol,0,bbrowij*sizeof(PetscInt));
    for(j1=0; j1<nreg; j1++) {
      MatGetColumnNorms(submatCij[j1+j3*(nreg+1)][0],NORM_1,ccolnorms);
      for(i=0; i<bbrowij; i++)if(ccolnorms[i]>0)bivinzcol[i]++;
      Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatBij[j1+j3*(nreg+1)][0]->data;//*aa=subA->data;
      ai= ab->i;
      nz=ab->nz;
      nrow=submatBij[j1+j3*(nreg+1)][0]->rmap->n;
      //printf("j1 %d nrow %d\n",j1,nrow);
      for(i=0; i<nrow-1; i++)if(ai[i]!=ai[i+1])bivinzrow[i]++;
      if(ai[i]<nz)bivinzrow[i]++;
    }
    ddrowi=submatBij[j3*(nreg+1)][0]->rmap->n;
    //printf("ddrowi %d\n",ddrowi);
    for(i=0; i<ddrowi; i++) {
      if(bivinzrow[i]>0) {
        bivinzrow[i]=bivirowsize;
        bivirowsize++;
      }
      if(bivinzcol[i]>0) {
        bivinzcol[i]=bivicolsize;
        bivicolsize++;
      }
    }
    //printf("j3 %d\n",j3);
    bivirowsize--;
    bivicolsize--;
    bivinzrow1=realloc (bivinzrow1,bivirowsize*sizeof(long int));
    bivinzcol1=realloc (bivinzcol1,bivicolsize*sizeof(int));
    j=0;
    j2=0;
    for(i=0; i<ddrowi; i++) {
      if(bivinzrow[i]>0) {
        bivinzrow1[j]=i*bbrowij;
        j++;
      }
      if(bivinzcol[i]>0) {
        bivinzcol1[j2]=i;
        j2++;
      }
      bivinzrow[i]=(bivinzrow[i]-1)*bivicolsize;//faster MSOL
    }
    vecbivisize=((long int)bivirowsize)*((long int)bivicolsize);
    printf("vecbivisize %ld bivirowsize %d bivicolsize %d\n",vecbivisize,bivirowsize,bivicolsize);
    if(vecbivisize<0)printf("Overflow detected, the interface problem is too big!!!!!!\n");
    //vecbivi=realloc (vecbivi,vecbivisize*sizeof(int));
    ha_cgetype *vecbivi= (ha_cgetype *) calloc (vecbivisize,sizeof(ha_cgetype));
    //printf("vecbivisize %ld rank %d\n",vecbivisize,rank);
    xi1indx=0;
    for(j1=0; j1<nreg; j1++) {
      j4=j1+j3*(nreg+1);
      time(&timestr);//=clock();
      Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatAij[j4]->data;//*aa=subA->data;
      MatTranspose(submatCij[j4][0],MAT_INITIAL_MATRIX,&submatCT);
      MatDestroy(&submatCij[j4][0]);//1
      Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatCT->data;//*aa=subA->data;
      ai= aa->i;
      aj= aa->j;
      vals=aa->a;
      nz=aa->nz;
      nrow=submatAij[j4]->rmap->n;
      ncol=submatAij[j4]->cmap->n;
      lasize=ceil((laA/100.0)*nz);
      int *irn=(int *) calloc (lasize,sizeof(int));
      int *jcn=(int *) calloc (lasize,sizeof(int));
      ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
      insize[j4*insizes+13]=bivirowsize;
      insize[j4*insizes+14]=bivicolsize;
      insize[j4*insizes]=nrow;
      insize[j4*insizes+1]=ncol;
      insize[j4*insizes+2]=nz;
      aic= ac->i;
      ajc= ac->j;
      valsc=ac->a;
      nzc=ac->nz;
      nrowc=submatCT->rmap->n;
      ncolc=submatCT->cmap->n;
      insize[j4*insizes+3]=nrowc;
      insize[j4*insizes+4]=ncolc;
      memcpy (irn+nz,ai,(nrow+1)*sizeof(PetscInt));
      memcpy (jcn,aj,nz*sizeof(PetscInt));
      memcpy (values,vals,nz*sizeof(ha_cgetype));
      //printf("nz %d j %d\n",nz,j);
      insize[j4*insizes+2]=nz;
      insize[j4*insizes+5]=nzc;
      MatDestroy(&submatAij[j4]);//1
      insize[j4*insizes+9]=laA;
      insize[j4*insizes+10]=rank;
      insize[j4*insizes+11]=j4;
      MatGetSize(submatBij[j4][0],&nrowb,&ncolb);
      Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatBij[j4][0]->data;//*aa=subA->data;
      ai= ab->i;//irnbs+j1irnbs;
      aj= ab->j;//jcnbs+j1nz;
      vals=ab->a;//valbs+j1nz;
      nz=ab->nz;//nzbs[j1];
      insize[j4*insizes+6]=nrowb;//ncolb;
      insize[j4*insizes+7]=ncolb;//nrowb;
      insize[j4*insizes+8]=nz;
      prep48_msol_(insize+j4*insizes,irn,jcn,values,aic,ajc,valsc,ai,aj,vals,vecbivi,bivinzrow,bivinzcol);
      MatDestroy(&submatCT);
      MatDestroy(&submatBij[j4][0]);//1
      //free(yi1[j1]);
      free(irn);//1
      free(jcn);//1
      free(values);//1
      time(&timeend);
      //printf("t %d Submatrix %d rank %d calculation time %f\n",j3,j1,rank,difftime(timeend,timestr));
    }

    //printf("HHHEERE\n");

    j4=nreg+j3*(nreg+1);
    Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatAij[j4]->data;//*aa=subA->data;
    ai= aa->i;
    aj= aa->j;
    vals=aa->a;
    nz=aa->nz;

//    vecbivisize=((long int)bivirowsize)*((long int)bivicolsize);//(VecSize-sumrowcolin)*(VecSize-sumrowcolin);  //
    j=0;
    for(i=vecbivisize-1; i!=-1; --i)if(vecbivi[i]!=0)j++;
    int *irn1=(int *) calloc (nz+j,sizeof(int));
    int *jcn1=(int *) calloc (nz+j,sizeof(int));

    //if(j==0)j=1;
//    obiviindx1=realloc(obiviindx1,(j+1)*sizeof(long int));
//    biviindx1=obiviindx1;
//    biviindx1++;
    j2=0;
    for(i=0; i<bivirowsize; i++) {
      bivinzrow1[i]/=bbrowij;
      bivinzrow1[i]++;
    }
    for(i=0; i<bivicolsize; i++) {
      bivinzcol1[i]++;
    }
    //printf("HHHEERE nz %d j %d vs %ld\n",nz,j,vecbivisize);
    for(i=0; i<vecbivisize; ++i) {
      if(vecbivi[i]!=0) {
        vecbivi[j2]=vecbivi[i];
        irn1[j2]=bivinzrow1[i/bivicolsize];
        jcn1[j2]=bivinzcol1[i%bivicolsize];
        j2++;
      }
    }
    //printf("HHHEERE\n");
    bivinzrow1=realloc (bivinzrow1,1*sizeof(long int));
    bivinzcol1=realloc (bivinzcol1,1*sizeof(int));
    //if(j2==0){j2=1;biviindx1[0]=0;}
//    vecbivi=realloc(vecbivi,j2*sizeof(ha_cgetype));
    //printf("J2 %d\n",j2);
    //for(i=0;i<j2;i++)printf("aii %f indx %d\n",vecbivi[i],biviindx1[i]);
    nz0=j2;
    nrow=submatAij[j4]->rmap->n;
    ncol=submatAij[j4]->cmap->n;
    if(vecbivisize<nz+j)vecbivi=realloc(vecbivi,(nz+j)*sizeof(ha_cgetype));
    printf("D ncol %d nrow %d nz %d\n",ncol,nrow,nz);
//    obiviindx0=realloc(obiviindx0,(nz+1)*sizeof(long int));
//    biviindx0=obiviindx0;
//    biviindx0++;
//    j2=0;
    for(i=1; i<nrow; i++) {
      for(j=ai[i-1]; j<ai[i]; j++) {
        vecbivi[j2]=vals[j];
        irn1[j2]=i;
        jcn1[j2]=aj[j]+1;
        j2++;
      }
    }
    for(j=ai[i]; j<nz; j++) {
      vecbivi[j2]=vals[j];
      irn1[j2]=i;
      jcn1[j2]=aj[j]+1;
      j2++;
    }
    //printf("in C %ld\n",biviindx0[1]);
    vecbivi=realloc(vecbivi,j2*sizeof(ha_cgetype));
    nz1=j2;
//    my_spar_compl_(biviindx1,&nz0,biviindx0,&nz,&nz1);
//    printf("rank %d nz0 %d nz %d nz1 %d\n",rank,nz0,nz,nz1);
//    vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
//    printf("nz1 %d\n",nz1);
//    //for(j=0;j<nz;j++)printf("rank %d j %d bivi %ld v %lf\n",rank,j,obiviindx0[j],vals[j]);
//    my_spar_add1l_(vecbivi,obiviindx1,irn1,jcn1,&nz0,vals,obiviindx0,&nz,&nz1,&ncol);
//    printf("nz1 %d\n",nz1);
    MatDestroy(&submatAij[j4]);//1
//    printf("Arank %d nz0 %d nz %d nz1 %d\n",rank,nz0,nz,nz1);
    //begin save for sol
    if(j3<10)strcpy(j1name,"000");
    if(j3<100&&j3>9)strcpy(j1name,"00");
    if(j3<1000&&j3>99)strcpy(j1name,"0");
    if(j3>=1000)j1name[0]='\0';
    sprintf(filename, "%d",j3);
    strcat(j1name,filename);
    strcpy(filename,"_bivi");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "wb"))==NULL) {
      printf("Cannot open file.\n");
    }
    fwrt=fwrite(vecbivi, sizeof(ha_cgetype),nz1, presolfile);
    if(fwrt==0)printf("File write error.");
    fclose(presolfile);
    strcpy(filename,"_rbvi");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "wb"))==NULL) {
      printf("Cannot open file.\n");
    }
    fwrt=fwrite(irn1, sizeof(int), nz1, presolfile);
    if(fwrt== 0) printf("File write error.");
    fclose(presolfile);
    strcpy(filename,"_cbvi");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "wb"))==NULL) {
      printf("Cannot open file.\n");
    }
    fwrt=fwrite(jcn1, sizeof(int), nz1, presolfile);
    //if(rank==2)for(i=0;i<ncol;i++)printf("rank %d j1 %d col %d\n",rank,j3,jcn1[i]);
    if(fwrt== 0) printf("File write error.");
    fclose(presolfile);
    //end save for sol

    nz=nz1;
    obiviindx1=realloc(obiviindx1,sizeof(long int));
    obiviindx0=realloc(obiviindx0,sizeof(long int));
    ldsize=ceil((laDi/100.0)*nz1);
    irn1=realloc(irn1,ldsize*sizeof(int));
    jcn1=realloc(jcn1,ldsize*sizeof(int));
    vecbivi=realloc(vecbivi,ldsize*sizeof(ha_cgetype));
    insize[j4*insizes]=nrow;
    insize[j4*insizes+1]=ncol;
    insize[j4*insizes+2]=nz;
    insize[j4*insizes+10]=rank;
    insize[j4*insizes+11]=j4;
    insize[j4*insizes+9]=laDi;
    //for(i=0;i<nz;i++)printf("aii %f irn1 %d jcn %d\n",vecbivi[i],irn1[i],jcn[i]);
    //printf("ok!!!!!!!! vecbiuisize %d nrow %d\n",ncol,nrow);
    if(presol) {
      int *insized=(int *) calloc (5,sizeof(int));
      insized[0]=nrow;
      insized[1]=ncol;
      insized[2]=nz;
      insized[4]=laDi;
      //if(rank==2)for(i=0;i<ncol;i++)printf("arank %d j1 %d col %d sct %d svb %d\n",rank,j3,jcn[i],(int)sizeof(PetscReal),(int)sizeof(ha_cgetype));
      if(cntl6==0&&SORD==0)cntl6in=0.3;
      else cntl6in=cntl6;
      spec51_rank_(insized,&cntl6in,irn1,jcn1,vecbivi);
      //insized[3]=(int)(0.9*insized[3]);
      insized=realloc(insized,(5+nreg*insizes)*sizeof(int));
      for(i=0; i<nreg; i++) {
        for(j=0; j<insizes; j++) {
          insized[i*insizes+5+j]=insize[(i+j3*(nreg+1))*insizes+j];
        }
      }
      if(j3<10)strcpy(j1name,"000");
      if(j3<100&&j3>9)strcpy(j1name,"00");
      if(j3<1000&&j3>99)strcpy(j1name,"0");
      if(j3>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j3);
      strcat(j1name,filename);
      strcpy(filename,"_rank");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      printf("Rank %d RRRRR! ncol %d nrow %d rank %d nz %d cntl6in %lf\n",rank,nrow,ncol,insized[3],nz,cntl6in);
      fwrt=fwrite(insized, sizeof(int), 5+nreg*insizes, presolfile);
      if(fwrt==0)printf("File write error.");
      fclose(presolfile);
      strcpy(filename,"_row");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fwrite(irn1, sizeof(int), nrow, presolfile);
      if(fwrt== 0) printf("File write error.");
      fclose(presolfile);
      strcpy(filename,"_col");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fwrite(jcn1, sizeof(int), ncol, presolfile);
      //if(rank==2)for(i=0;i<ncol;i++)printf("rank %d j1 %d col %d\n",rank,j3,jcn1[i]);
      if(fwrt== 0) printf("File write error.");
      fclose(presolfile);
      free(insized);//1
    }
    else {
      prep48_alu_(insize+j4*insizes,irn1,jcn1,vecbivi);
    }
    //printf("ok!!!!!!!!\n");
    //vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
    free(vecbivi);//1
    free(irn1);//1
    free(jcn1);//1
  }
  //if(presol)return 0;
  free(ccolnorms);//1
  ierr = PetscFree(submatAij);//1
  printf("Rank %d End of preparation!\n",rank);
  MPI_Barrier(PETSC_COMM_WORLD);

  //bbrowij=submatB[0]->rmap->n;
  //printf("j3 %d bbrowij %d\n",j3,bbrowij);
  //ha_cgetype *vecbivi= (ha_cgetype *) calloc (vecbivisize,sizeof(ha_cgetype));
//  for(j1=0; j1<nmatint; j1++) {
//    MatDestroy(&submatCT);
//    free(yi1[j1]);//1
//    MatDestroy(&submatB1[j1][0]);//submatBT);//1
//    ierr = PetscFree(submatB1[j1]);//1
//    MatDestroy(&submatB2[j1][0]);//submatBT);//1
//    ierr = PetscFree(submatB2[j1]);//1
  //ierr = PetscFree(ai12);
  //ierr = PetscFree(aj12);
  //ierr = PetscFree(vals12);
//  }
  free(bivinzcol);//1
  free(bivinzrow);//1
//  ierr = PetscFree(submatB1);//1
//  ierr = PetscFree(submatB2);//1
//  CHKERRQ(ierr);
//  free(yi1);//1
  free(bivinzrow1);//1
  free(bivinzcol1);//1
  //vecbivi=realloc(vecbivi,j2*sizeof(ha_cgetype));
//  if(rank==mpisize-1) {
//    MatDestroy(&submatD[0]);//1
//    //free(irn1);
//    //free(jcn);
//    //free(insizeD);
//    //x0=(ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));//realloc (x0,*sizeof(ha_cgetype));
//  }
//  free(vecbiui);//1
  free(obiviindx1);//1
  free(obiviindx0);//1
  //free(vecbivi);//1
//  free(indicesC);//1
//  PetscFree(submatD);//1
//  for(j1=0; j1<nmatint; j1++) {
//    for(i=0; i<nreg; i++) {
//      MatDestroy(&submatCij[j1*(nreg+1)+i][0]);//1
//      MatDestroy(&submatBij[j1*(nreg+1)+i][0]);//1
//    }
//  }
  for(i=0; i<nmatin; i++) {
    PetscFree(submatCij[i]);//1
    PetscFree(submatBij[i]);//1
  }
//  ierr = PetscFree(colindicesbc1);//1
//  ierr = PetscFree(colindicesbc2);//1
  PetscFree(submatCij);//1
  PetscFree(submatBij);//1

//  free(xi1);//1
  free(insize);//1
  //ierr = PetscFree(submatC);//1
  free(begblock);//1
  timestr=clock();
  for(j1=0; j1<nmatin; j1++) {
    free(fn01[j1]);//1
    free(fn02[j1]);//1
    free(fn03[j1]);//1
  }
  free(fn01);//1
  free(fn02);//1
  free(fn03);//1
  free(cirank);//1
  free(birank);//1
  free(aord);//1
  free(ndblockinrank);//1
  return 0;
}

int HaNDBBDParSol(Mat A, Vec b, ha_cgetype *x1, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, uvadd ndblock,uvadd nreg,uvadd ntime,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,uvdim laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol) {//,bool iter
  IS *rowindices=NULL,*colindices=NULL,*rowindicesbc=NULL,*colindicesbc1=NULL,*colindicesbc2=NULL,*Cindices=NULL,*Bindices=NULL,Cindicesc,Bindicesc,*BBindices=NULL,*BBiindices=NULL,*rowBBij=NULL,*colBBij=NULL;//,*colindicesbcpm,*colindicesbcpm1
  const PetscInt *nindices;
  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nmatint,nmatinplust,nmatminust,nrowcolin,sumrowcolin,i,i1,j,j0,j1,j2,j3,j4,j5,j6,j7,l0,l1,l2,l3,l4,l5,rank1,proc1=0,nnzmax,j1nz,j1irnbs;
  Mat *submatAij=NULL,*submatC=NULL,*submatD=NULL,*submatBB=NULL,*submatBBij=NULL;//,*submatCij,*submatBij;,*submatB
  Vec vecxd;
  PetscInt *ai,*aj,*a1i,*a1j,*a2i,*a2j,*aic,*ajc;
  PetscReal cntl6in;
  uvadd ldsize;
  long int *biviindx1,*biviindx0,j7l;
  long int *obiviindx1= NULL;//(long int*)calloc(1,sizeof(long int));
  long int *obiviindx0= NULL;//(long int*)calloc(1,sizeof(long int));
  size_t freadresult,fwrt,frd;
  //printf("l l i %d l i %d i %d",sizeof(long long int),sizeof(long int),sizeof(int));
//  char tempchar[256];
  //obiviindx1=realloc(obiviindx1,sizeof(int));
  //obiviindx1[0]=-1;
  //obiviindx0=realloc(obiviindx0,sizeof(int));
  //obiviindx0[0]=-1;
  PetscInt nrow,ncol,nz,nrowc,ncolc,nrowb,ncolb,nzc,nzv;
  //forint iv,jv,la1;
  MPI_Status   status;
  clock_t timestr,timeend;
  //forint forintwrp;
  int j2int,la1,colcut,rowcut,jthrd,nthrd=1,unequal=0;
  nthrd=omp_get_max_threads( );
  bool ifremove=false,iterstop;
  //long long int nrow,ncol,nz,nrowc,ncolc,nzc,nzv;
  PetscScalar *vals,*val1s,*val2s,*valsc,vecval;
  ha_cgetype vval;
  PetscErrorCode ierr;
  PetscViewer viewer;
  FILE* fp1,*fp2,*fp3,*presolfile;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  int *ndblockinrank= (int *) calloc (mpisize,sizeof(int));
  char filename[1024],rankname[1024],j1name[1024],tempchar[1024];//,fn1[1024],fn2[1024],fn3[1024]
  if(rank<10)strcpy(rankname,"000");
  if(rank<100&&rank>9)strcpy(rankname,"00");
  if(rank<1000&&rank>99)strcpy(rankname,"0");
  if(rank>=1000)rankname[0]='\0';
  sprintf(filename, "%d",rank);
  strcat(rankname,filename);
  int nmatinBB=1;
  //reg matrix in one process only
  //nmatin=(uvadd)ndblock/mpisize;
  //if(nmatin*mpisize<ndblock)nmatinplus=nmatin+1;
  //else nmatinplus=nmatin;
  nmatint=(uvadd)ntime/mpisize;
  if(nmatint*mpisize<ntime)unequal=1;
  nmatminust=nmatint;
  if(nmatint*mpisize<ntime)nmatinplust=nmatint+1;
  else nmatinplust=nmatint;
  for(i=0; i<mpisize; i++)if(rank+1<=ntime-mpisize*nmatint)nmatint++;
  nmatin=(nreg+1)*nmatint;
  nmatinplus=(nreg+1)*nmatinplust;
  begblock[rank]=nmatin;
  printf("rank %d nmatin %d nmatint %d nmplus %d\n",rank,nmatin,nmatint,nmatinplus);
  for(i=0; i<mpisize; i++) {
    j=rank;
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    j=j2;
    MPI_Bcast(&begblock[j2],1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  j2=begblock[0];
  begblock[0]=0;
  for(i=1; i<mpisize; i++) {
    j1=begblock[i];
    begblock[i]=begblock[i-1]+j2;
    j2=j1;
  }
  for(i=0; i<mpisize; i++)printf("rank %d beg block %d\n",rank,begblock[i]);
  j2=nmatin;
  MPI_Bcast(&j2,1, MPI_INT,0, PETSC_COMM_WORLD);
  ndblockinrank[0]=j2;
  for(i=1; i<mpisize; i++) {
    j2=nmatin;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    ndblockinrank[i]=ndblockinrank[i-1]+j2;
  }
  //int *birank= (int *) calloc (nmatin,sizeof(int));
  //int *cirank= (int *) calloc (nmatin,sizeof(int));
  //for(i=0; i<nmatin; i++)cirank[i]=-1;
  //int *aord= (int *) calloc (nmatin,sizeof(int));
  for(i=0; i<ndblock; i++) {
    j2=i-ndblockinrank[rank];
    if(j2<0&&j2>=-nmatin) {
      j1=i/(nreg+1);
      j=i-(nreg+1)*j1;
      if(rank==0)j3=0;
      else j3=ndblockinrank[rank-1];
      //if(j<nreg)cirank[j2+nmatin]=(nreg+1)*j1+nreg-j3;
      //else cirank[j2+nmatin]=-1;
    }
  }
  //printf("OK1!!!\n");
  j1=0;
  for(j=0; j<mpisize; j++)if(j<rank)j1+=ndblockinrank[j];
  free(ndblockinrank);
//   for(i=0; i<nmatin; i++) {
//     j2=(j1+i)/(nreg+1);
//     if(i-j2*(nreg+1)!=0)aord[i]=1;
//   }
  //printf("OK1a!!!\n");
  ierr = PetscMalloc(nmatin*sizeof(IS **),&rowindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&rowindicesbc);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&colindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbc1);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbc2);
  CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbcpm);
  //CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbcpm1);
  //CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&Cindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&Bindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatinBB*sizeof(IS **),&BBindices);
  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatAij);
//  CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatint*sizeof(Mat *),&submatC);
  //CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatCij);
  //CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatint*sizeof(Mat *),&submatB);
  //CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatBij);
  //CHKERRQ(ierr);
  //printf("OK2!!!\n");
  //ierr = PetscMalloc(nmatinBB*sizeof(Mat *),&submatBB);
  //CHKERRQ(ierr);
  Mat** submatBij = (Mat**)calloc(nmatin,sizeof(Mat*));
//   for (i=0; i<nmatin; i++) {
//     submatBij[i] = (Mat*)calloc(1,sizeof(Mat));
//   }
  Mat** submatCij = (Mat**)calloc(nmatin,sizeof(Mat*));
//   for (i=0; i<nmatin; i++) {
//     submatCij[i] = (Mat*)calloc(1,sizeof(Mat));
//   }
  Mat** submatB1 = (Mat**)calloc(nmatint,sizeof(Mat*));
  for (i=0; i<nmatint; i++) {
    submatB1[i] = (Mat*)calloc(1,sizeof(Mat));
  }
  Mat** submatB2 = (Mat**)calloc(nmatint,sizeof(Mat*));
  for (i=0; i<nmatint; i++) {
    submatB2[i] = (Mat*)calloc(1,sizeof(Mat));
  }
  //printf("rank %d OK3!!!\n",rank);

  ha_cgetype **yi2= (ha_cgetype**)calloc(nmatin,sizeof(ha_cgetype*));
  for (i=0; i<nmatin; i++) {
    //printf("i %d i+b %d s %d\n",i,i+begblock[rank],ha_ndblocks[i+begblock[rank]]);
    yi2[i] = (ha_cgetype*)calloc(ha_ndblocks[i+begblock[rank]],sizeof(ha_cgetype));
  }
  //printf("rank %d OK4!!!\n",rank);
  sumrowcolin=0;
  for(j=0; j<ndblock; j++) {
    sumrowcolin+=ha_ndblocks[j];
    //printf("wwwwwww myid %d j %d block %d\n",rank,j,ha_ndblocks[j]);
  }
  printf("Vecs %ld sumrow %d!!!\n",VecSize,sumrowcolin);
  //PetscScalar *yd = (PetscScalar*)calloc(VecSize-sumrowcolin,sizeof(PetscScalar));
  PetscInt *offblock= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
  PetscInt *offblockrow= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
  //printf("rank %d OK6!!!\n",rank);
  offblock[0]=0;
  for(j=0; j<ndblock; j++)offblock[j+1]=offblock[j]+countvarintra1[j+1]-countvarintra1[j]-ha_ndblocks[j];
  //printf("rank %d OK7!!!\n",rank);
  offblockrow[0]=0;
  for(j=0; j<ndblock; j++)offblockrow[j+1]=offblockrow[j]+counteqnoadd[j]-ha_ndblocks[j];
  //printf("OK123\n");
  //printf("rank %d OK8!!!\n",rank);

  char **fn01= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn01[i] = (char*)calloc(1024,sizeof(char));
  char **fn02= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn02[i] = (char*)calloc(1024,sizeof(char));
  char **fn03= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn03[i] = (char*)calloc(1024,sizeof(char));
  //printf("rank %d OK9!!!\n",rank);

  timestr=clock();
  //printf("rank %d OK10!\n",rank);
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      if(j1<10)strcpy(j1name,"000");
      if(j1<100&&j1>9)strcpy(j1name,"00");
      if(j1<1000&&j1>99)strcpy(j1name,"0");
      if(j1>=1000)j1name[0]='\0';
      //printf("rank %d nmatin %d matinp %d j1 %d\n",rank,nmatin,nmatinplus,j1);
      sprintf(filename, "%d",j1);
      strcat(j1name,filename);
      strcpy(filename,"_irnv");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn01[j1],filename);
      strcpy(filename,"_keep");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn02[j1],filename);
      strcpy(filename,"_vav");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn03[j1],filename);
    }
  }
  printf("rank1w2 %d\n",rank);
  //printf("rank %d OK11\n",rank);
  //int longi=3*pow(10,9);
  //printf("OK012345 rank %d mpis %d sizeof int %ld size of long int %ld\n",rank,mpisize,sizeof(int),sizeof(long int));
  for(j=0; j<mpisize; j++) {
    //printf("j %d rank %d\n",j,rank);
    for(j3=0; j3<nmatinplus; j3++) {
      if(j3>=nmatin)i=nmatin-1;
      else i=j3;
      //printf("rank %d OK1!!! i %d j %d nmatin %d\n",rank,i,j,nmatin);
      //printf("rank1wc %d i %d j3 %d i %d nmatin %d beg %d\n",rank,i+begblock[j],j3,i,nmatin,begblock[j]);
      bfirst=counteq[i+begblock[rank]];
      bend=ha_ndblocks[i+begblock[rank]];
      //printf("rank1wc %d i %d\n",rank,i+begblock[j]);
      MPI_Bcast(&bfirst,1, MPI_INT,j, PETSC_COMM_WORLD);
      MPI_Bcast(&bend,1, MPI_INT,j, PETSC_COMM_WORLD);
      PetscInt *indices= (PetscInt *) calloc (bend,sizeof(PetscInt));
      ha_cgetype *yi0 = (ha_cgetype*)calloc(bend,sizeof(ha_cgetype));
      for(j1=0; j1<bend; j1++) { //ha_ndblocks[i+begblock[j]]
        j2=ha_rows[bfirst+j1];//+bfirst;
        if(j2>=Istart&&j2<Iend)indices[j1]=j2;
        else indices[j1]=-1;
      }
      //printf("rank %d i %d j %d bf %d be %d\n",rank,i,j,bfirst,bend);
      VecGetValues(b,bend,indices,yi0);
      HaReduce(yi0,(forint)bend,mpisize,rank,j);
      if(rank==j&&i<nmatin) {
        memcpy (yi2[i],yi0,bend*sizeof(ha_cgetype));
      }
      free(indices);
      free(yi0);
    }
  }
  printf("rank1w2 %d\n",rank);
  ha_cgetype **yi1= (ha_cgetype**)calloc(nmatint,sizeof(ha_cgetype*));
  for (i=0; i<nmatint; i++) {
    j1=0;
    for (j=0; j<nreg+1; j++) {
      j1+=ha_ndblocks[i*(nreg+1)+j+begblock[rank]];
    }
    yi1[i] = (ha_cgetype*)calloc(j1,sizeof(ha_cgetype));
  }
  for (i=0; i<nmatint; i++) {
    j1=0;
    for (j=0; j<nreg+1; j++) {
      memcpy (&yi1[i][j1],yi2[i*(nreg+1)+j],ha_ndblocks[i*(nreg+1)+j+begblock[rank]]*sizeof(ha_cgetype));
      free(yi2[i*(nreg+1)+j]);
      j1+=ha_ndblocks[i*(nreg+1)+j+begblock[rank]];
    }
  }
  free(yi2);
  //printf("OK12345 rank %d\n",rank);
  j1=0;
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (ha_ndblocks[i+begblock[rank]],sizeof(PetscInt));
    //printf("ha_row 0 %d 588 %d\n",ha_rows[0],ha_rows[588]);
    //if(rank==0)printf("rank %d i %d bfirst %d block %d\n",rank,i,bfirst,ha_ndblocks[i+begblock[rank]]);
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_rows[bfirst+j];//+bfirst;
      //if(rank==0&&i==0)printf("indices %d rank %d\n",indices[j],ha_ndblocks[i+begblock[rank]]);
    }

//    bend=counteq[i+begblock[rank]]+counteqnoadd[i+begblock[rank]];
//    for(j=offblockrow[i+begblock[rank]]; j<offblockrow[i+1+begblock[rank]]; j++) {
//      j1=j-offblockrow[i+begblock[rank]];
//      indicesB[j]=ha_rows[bend-j1-1]+bfirst;
//    }
    //if(rank==0)printf("rank %d i %d ha_ndblocks[i+begblock[rank]] %d\n",rank,i,ha_ndblocks[i+begblock[rank]]);
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
    //if(cirank[i]>-1)ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
    bfirst=countvarintra1[i+begblock[rank]];
    //if(rank==0)printf("rank %d i %d bfirst %d block %d\n",rank,i,bfirst,ha_ndblocks[i+begblock[rank]]);
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_cols[bfirst+j];//+bfirst;
      //if(rank==0&&i==90)printf("indices %d rank %d\n",indices[j],ha_ndblocks[i+begblock[rank]]);
    }
//    bend=countvarintra1[i+1+begblock[rank]];
//    for(j=offblock[i+begblock[rank]]; j<offblock[i+1+begblock[rank]]; j++) {
//      j1=j-offblock[i+begblock[rank]];
//      indicesC[j]=ha_cols[bend-j1-1]+bfirst;
//      if(rank==1)printf("rank %d cindices %d bend %d j1 %d j %d\n",rank,indicesC[j],ha_cols[bend-j1-1],j1,j);
//    }
    //MPI_Bcast(indicesC[offblock[i+begblock[rank]]],(offblock[i+1+begblock[rank]]-offblock[i+begblock[rank]])*sizeof(PetscInt), MPI_BYTE,rank, PETSC_COMM_WORLD);
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,colindices+i);
    //if(cirank[i]>-1)ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,colindicesij+i);
    free(indices);
  }
  //MPI_Allreduce(MPI_IN_PLACE,indicesC,VecSize-sumrowcolin, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
  for(i=0; i<nmatint; i++) {
    j2=0;
    for(j1=0; j1<nreg+1; j1++)j2+=ha_ndblocks[i*(nreg+1)+j1+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (j2,sizeof(PetscInt));
    j2=0;
    //j3=0;
    for(j1=0; j1<nreg+1; j1++) {
      bfirst=counteq[i*(nreg+1)+j1+begblock[rank]];
      //printf("rank %d bfirst %d\n",rank,bfirst);
      //printf("ha_row 0 %d 588 %d\n",ha_rows[0],ha_rows[588]);
      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
        indices[j2]=ha_rows[bfirst+j];//+bfirst;
        //if(rank==0&&i==nmatint-1)printf("rank %d indx %d\n",rank,indices[j2]);
        j2++;
        //if(rank==0&&i==1)printf("indices %d rank %d\n",indices[j],ha_ndblocks[i+rank*nmatinplus]);
      }
    }
    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,rowindicesbc+i);
    //if(cirank[i]>-1)ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
//    j2=0;
//    for(j1=0; j1<nreg+1; j1++) {
//      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
//      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
//        indices[j2]=j2+bfirst;
//        j2++;
//      }
//    }
//    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbcpm+i);
//    j2=0;
//    for(j1=0; j1<nreg+1; j1++) {
//      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
//      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
//        indices[j2]=ha_cols[bfirst+j]-bfirst;//+bfirst;
//        //if(i==0&&j1==2)printf("j2 %d cindices %d cindx-1 %d\n",j2,indices[j2],indices[j2-1]);
//        //if(j2>0)if(indices[j2]<=indices[j2-1])printf("i %d j1 %d j2 %d indx %d\n",i,j1,j2,indices[j2]);
//        j2++;
//      }
//    }
//    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbcpm1+i);
    j2=0;
    for(j1=0; j1<nreg; j1++) {
      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
        indices[j2]=ha_cols[bfirst+j];//+bfirst;
        //if(i==0&&j1==2)printf("j2 %d cindices %d cindx-1 %d\n",j2,indices[j2],indices[j2-1]);
        //if(j2>0)if(indices[j2]<=indices[j2-1])printf("i %d j1 %d j2 %d indx %d\n",i,j1,j2,indices[j2]);
        j2++;
      }
    }
    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbc1+i);
    j2=0;
    for(j1=nreg; j1<nreg+1; j1++) {
      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
        indices[j2]=ha_cols[bfirst+j];//+bfirst;
        //if(i==0&&j1==2)printf("j2 %d cindices %d cindx-1 %d\n",j2,indices[j2],indices[j2-1]);
        //if(j2>0)if(indices[j2]<=indices[j2-1])printf("i %d j1 %d j2 %d indx %d\n",i,j1,j2,indices[j2]);
        j2++;
      }
    }
    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbc2+i);
    free(indices);
  }
  //printf("OK12345 rank %d\n",rank);
  //ISCopy(colindices,colindicesij);
  //ISCopy(rowindices,rowindicesij);
  //printf("0okkkk\n");
  PetscInt *indicesC= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
  PetscInt *indicesB= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
  for(i=0; i<ndblock; i++) {
    bend=countvarintra1[i+1];
    //bfirst=countvarintra1[i];
    for(j=offblock[i]; j<offblock[i+1]; j++) {
      j1=j-offblock[i];
      indicesC[j]=ha_cols[bend-j1-1];//+bfirst;
      //if(i==98)printf("rank %d cindices %d bend %d j1 %d j %d\n",rank,indicesC[j],ha_cols[bend-j1-1],j1,j);
    }
  }
  j1=offblock[ndblock];
  //printf("rank %d j1 %d count %d vecs %d ndblock %d v-s %d\n",rank,j1,countvarintra1[ndblock],VecSize,ndblock,VecSize-sumrowcolin);
  for(i=countvarintra1[ndblock]; i<VecSize; i++) {
    indicesC[j1]=i;
    //printf("rank %d i %d indx %d\n",rank,i,indicesC[j1]);
    j1++;
  }
  //printf("1okkkk\n");
  for(i=0; i<ndblock; i++) {
    //bfirst=counteq[i];
    bend=counteq[i]+counteqnoadd[i];
    for(j=offblockrow[i]; j<offblockrow[i+1]; j++) {
      j1=j-offblockrow[i];
      indicesB[j]=ha_rows[bend-j1-1];//+bfirst;
    }
  }
  j1=offblockrow[ndblock];
  for(i=counteq[ndblock-1]+counteqnoadd[ndblock-1]; i<VecSize; i++) {
    indicesB[j1]=i;
    j1++;
  }
  //printf("2okkkk\n");
//  for(i=0; i<mpisize; i++) {
//    j=offblockrow[begblock[rank]];
//    j1=offblockrow[nmatin+begblock[rank]]-offblockrow[begblock[rank]];
//    MPI_Bcast(&j,1, MPI_INT,i, PETSC_COMM_WORLD);
//    MPI_Bcast(&j1,1, MPI_INT,i, PETSC_COMM_WORLD);
//    MPI_Bcast(&indicesB[j],j1, MPI_INT,i, PETSC_COMM_WORLD);
//  }
//  for(i=0; i<VecSize-sumrowcolin; i++)printf("rank %d indxb %d\n",rank,indicesB[i]);
  for(i=0; i<VecSize-sumrowcolin; i++)if(indicesB[i]>=Istart) {
      j=i;
      break;
    }
  for(j2=i; j2<VecSize-sumrowcolin; j2++)if(indicesB[j2]>=Iend) {
      j1=j2;
      break;
    }
  if(j2==VecSize-sumrowcolin)j1=VecSize-sumrowcolin;
  printf("rank %d j %d j1 %d istart %d iend %d\n",rank,j,j1,Istart,Iend);
  ha_cgetype *vecbiui= (ha_cgetype *) calloc (VecSize-sumrowcolin,sizeof(ha_cgetype));
  VecGetValues(b,j1-j,&indicesB[j],&vecbiui[j]);//implicite yd save mem
  ierr = VecDestroy(&b);
  CHKERRQ(ierr);
  //printf("rank2 %d\n",rank);
//  for(i=0; i<mpisize; i++) {
//    j=offblock[begblock[rank]];
//    j1=offblock[nmatin+begblock[rank]]-offblock[begblock[rank]];
//    MPI_Bcast(&j,1, MPI_INT,i, PETSC_COMM_WORLD);
//    MPI_Bcast(&j1,1, MPI_INT,i, PETSC_COMM_WORLD);
//    MPI_Bcast(&indicesC[j],j1, MPI_INT,i, PETSC_COMM_WORLD);
//  }
//  printf("OK132\n");
  ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesC,PETSC_COPY_VALUES,&Cindicesc);
  ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesB,PETSC_COPY_VALUES,&Bindicesc);
  for(i=0; i<nmatint; i++) {
    //ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesC,PETSC_COPY_VALUES,Cindices+i);
    //ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesB,PETSC_COPY_VALUES,Bindices+i);
    Cindices[i]=Cindicesc;
    Bindices[i]=Bindicesc;
  }
  free(indicesB);
  free(offblock);
  free(offblockrow);
  PetscInt *indexBB= (PetscInt *) calloc (VecSize,sizeof(PetscInt));
  for(i=0; i<VecSize; i++)indexBB[i]=i;
  ISCreateGeneral(PETSC_COMM_SELF,VecSize,indexBB,PETSC_COPY_VALUES,BBindices);
  printf("Begin preparation rank %d\n",rank);
//  ierr = MatCreateSubMatrices(A,nmatin,rowindices,colindices,MAT_INITIAL_MATRIX,&submatAij);
  //if(rank==0)printf("rank %d i %d bbrowij %d\n",rank,(nmatint-1)*(nreg+1)+nreg,submatAij[(nmatint-1)*(nreg+1)+nreg]->cmap->n);
  int insizes=17;
  int *insize=(int *) calloc (insizes*nmatin,sizeof(int));
//  for(j1=0; j1<nmatin; j1++){
//    if(aord[j1]==1){
//      time(&timestr);//=clock();
//      Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatA[j1]->data;//*aa=subA->data;
//      ai= aa->i;
//      aj= aa->j;
//      vals=aa->a;
//      nz=aa->nz;
//      nrow=submatA[j1]->rmap->n;
//      ncol=submatA[j1]->cmap->n;
//      int *irn=(int *) calloc (laA*nz,sizeof(int));
//      int *jcn=(int *) calloc (laA*nz,sizeof(int));
//      ha_cgetype *values= (ha_cgetype *) calloc (laA*nz,sizeof(ha_cgetype));
//      insize[j1*insizes]=nrow;
//      insize[j1*insizes+1]=ncol;
//      insize[j1*insizes+2]=nz;
//      memcpy (irn+nz,ai,(nrow+1)*sizeof(PetscInt));
//      memcpy (jcn,aj,nz*sizeof(PetscInt));
//      memcpy (values,vals,nz*sizeof(ha_cgetype));
//      MatDestroy(&submatA[j1]);
//      insize[j1*insizes+9]=laA;
//      insize[j1*insizes+10]=rank;
//      insize[j1*insizes+11]=j1;
//      prep48_alu_(insize+j1*insizes,irn,jcn,values);
//      free(irn);
//      free(jcn);
//      free(values);
//      time(&timeend);
//      printf("Submatrix %d rank %d calculation time %f\n",j1,rank,difftime(timeend,timestr));
//    }
//  }
  printf("rank %d A\n",rank);
  PetscInt *indicesbbij=(int *) calloc (1,sizeof(PetscInt));
  ierr = PetscMalloc(1*sizeof(IS **),&rowBBij);
  CHKERRQ(ierr);
  ierr = PetscMalloc(1*sizeof(IS **),&colBBij);
  CHKERRQ(ierr);

  for(i=0; i<nmatminust; i++) {
    j2=i*(nreg+1)+nreg;
    ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindicesbc1+i,MAT_INITIAL_MATRIX,&submatBBij);
    nrow=submatBBij[0]->rmap->n;
    ncol=submatBBij[0]->cmap->n;
    if(nrow>ncol) {
      indicesbbij=realloc(indicesbbij,nrow*sizeof(PetscInt));
      for(j=0; j<nrow; j++)indicesbbij[j]=j;
    }
    else {
      indicesbbij=realloc(indicesbbij,ncol*sizeof(PetscInt));
      for(j=0; j<ncol; j++)indicesbbij[j]=j;
    }
    ISCreateGeneral(PETSC_COMM_SELF,nrow,indicesbbij,PETSC_COPY_VALUES,rowBBij);
    j3=0;
    for(j=0; j<nreg; j++) {
      j1=i*(nreg+1)+j;
      ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
      ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[j1+begblock[rank]],indicesbbij+j3,PETSC_COPY_VALUES,colBBij);
      ierr = MatCreateSubMatrices(submatBBij[0],1,rowBBij,colBBij,MAT_INITIAL_MATRIX,&submatBij[j1]);
      ierr = ISDestroy(&colBBij[0]);
      j3+=ha_ndblocks[j1+begblock[rank]];
      //printf("rank %d i %d j3 %d bbrowij %d\n",rank,i,j,submatCij[j1][0]->cmap->n);
    }
    ierr = MatDestroy(&submatBBij[0]);
    ierr = PetscFree(submatBBij);//1
    submatBBij=NULL;
    ierr = ISDestroy(&rowBBij[0]);
  }
  //indicesbbij=realloc(indicesbbij,sizeof(PetscInt));
  //printf("rank %d A1\n",rank);
//  if(nmatin==nmatinplus)strcpy(tempchar,"MAT_INITIAL_MATRIX");
//  else strcpy(tempchar,"MAT_REUSE_MATRIX");
  if(unequal==1){
  if(nmatint<nmatinplust) for(i=0; i<nreg; i++) {//to create again
      MatDestroy(&submatCij[(nmatint-1)*(nreg+1)+i][0]);
      ierr = PetscFree(submatCij[(nmatint-1)*(nreg+1)+i]);
      MatDestroy(&submatBij[(nmatint-1)*(nreg+1)+i][0]);
      ierr = PetscFree(submatBij[(nmatint-1)*(nreg+1)+i]);
    }
  j2=(nmatint-1)*(nreg+1)+nreg;
  ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindicesbc1+(nmatint-1),MAT_INITIAL_MATRIX,&submatBBij);
  nrow=submatBBij[0]->rmap->n;
  ncol=submatBBij[0]->cmap->n;
  if(nrow>ncol) {
    indicesbbij=realloc(indicesbbij,nrow*sizeof(PetscInt));
    for(j=0; j<nrow; j++)indicesbbij[j]=j;
  }
  else {
    indicesbbij=realloc(indicesbbij,ncol*sizeof(PetscInt));
    for(j=0; j<ncol; j++)indicesbbij[j]=j;
  }
  ISCreateGeneral(PETSC_COMM_SELF,nrow,indicesbbij,PETSC_COPY_VALUES,rowBBij);
  j3=0;
  for(j=0; j<nreg; j++) {
    j1=(nmatint-1)*(nreg+1)+j;
    ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[j1+begblock[rank]],indicesbbij+j3,PETSC_COPY_VALUES,colBBij);
    ierr = MatCreateSubMatrices(submatBBij[0],1,rowBBij,colBBij,MAT_INITIAL_MATRIX,&submatBij[j1]);
    ierr = ISDestroy(&colBBij[0]);
    j3+=ha_ndblocks[j1+begblock[rank]];
    //printf("rank %d i %d j3 %d bbrowij %d\n",rank,i,j,submatCij[j1][0]->cmap->n);
  }
  ierr = MatDestroy(&submatBBij[0]);
  ierr = PetscFree(submatBBij);
  submatBBij=NULL;
  ierr = ISDestroy(&rowBBij[0]);
  }
  free(begblock);
  free(indicesbbij);
  PetscFree(rowBBij);
  PetscFree(colBBij);

//  for(j=0; j<nreg; j++) {
//    j1=(nmatint-1)*(nreg+1)+j;
//    j2=(nmatint-1)*(nreg+1)+nreg;
//    ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
//    ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindices+j1,MAT_INITIAL_MATRIX,&submatBij[j1]);
//    //printf("rank %d i %d j3 %d bbrowij %d\n",rank,i,j1,submatCij[j1][0]->cmap->n);
//  }
  printf("rank %d A\n",rank);
  ierr = MatCreateSubMatrices(A,nmatint,rowindicesbc,Cindices,MAT_INITIAL_MATRIX,&submatC);
  printf("C\n");
  ierr = MatCreateSubMatrices(A,nmatinBB,Bindices,BBindices,MAT_INITIAL_MATRIX,&submatBB);
  ierr = MatDestroy(&A);
  CHKERRQ(ierr);
  PetscInt BBrow,BBcol;//,CCrow,CCcol;
  MatGetSize(submatBB[0],&BBrow,&BBcol);
  //MatGetSize(submatC[0],&CCrow,&CCcol);
  printf("rank %d bbrow %d bbcol %d\n",rank,BBrow,BBcol);//,CCrow,CCcol);
  PetscInt *indexBBi= (PetscInt *) calloc (BBrow,sizeof(PetscInt));
  for(i=0; i<BBrow; i++)indexBBi[i]=i;
  ierr = ISDestroy(&Bindicesc);
  ISCreateGeneral(PETSC_COMM_SELF,BBrow,indexBBi,PETSC_COPY_VALUES,&Bindicesc);
  for(i=0; i<nmatint; i++) {
//    ierr = ISDestroy(&Bindices[i]);
//    ISCreateGeneral(PETSC_COMM_SELF,BBrow,indexBBi,PETSC_COPY_VALUES,Bindices+i);
    Bindices[i]=Bindicesc;
  }
  //ierr = MatCreateSubMatrices(submatBB[0],nmatint,Bindices,colindicesbc,MAT_INITIAL_MATRIX,&submatB);
  //char           ordering[256] = MATORDERINGRCM;
//         #pragma omp parallel private(i,ierr)
//         {
//         #pragma omp for
  for(i=0; i<nmatint; i++) {
    ierr = MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc1[i],MAT_INITIAL_MATRIX,&submatB1[i][0]);
    ierr = MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc2[i],MAT_INITIAL_MATRIX,&submatB2[i][0]);
    //ierr = MatCreateSubMatrices(submatBB[0],1,Bindices+i,colindicesbc+i,MAT_INITIAL_MATRIX,&submatB[i]);
    //ierr = MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc[i],MAT_INITIAL_MATRIX,&submatB[i][0]);
    //MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbcpm[i],MAT_INITIAL_MATRIX,&submatBT);
    //MatGetOrdering(A,ordering,&Bindices[i],&colindicesbcpm1[i]);
    //ISSetPermutation(Bindices[i]);
    //ISSetPermutation(colindicesbcpm1[i]);
    //MatPermute(submatBT,Bindices[i],colindicesbcpm1[i],&submatB[i][0]);
    //MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc[i],&submatBT);
    //MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc[i],&submatB[i][0]);
    //MatGetLocalSubMatrix(submatBB[0],Bindices[i],colindicesbc[i],&submatBT);
    //MatGetLocalSubMatrix(submatBB[0],Bindices[i],colindicesbc[i],&submatB[i][0]);
    //MatCopy(submatBT,submatB[i][0],DIFFERENT_NONZERO_PATTERN);
    //MatConvert(submatBT,MATSEQAIJ,MAT_INITIAL_MATRIX,&submatB[i][0]);
    //MatDuplicate(submatBT,MAT_COPY_VALUES,&submatB[i][0]);
    //MatRestoreLocalSubMatrix(submatBB[0],Bindices[i],colindicesbc[i],&submatBT);
    //printf("iiiii %d\n",i);
    //ierr = ISDestroy(&colindicesbcpm[i]);
    //ierr = ISDestroy(&colindicesbcpm1[i]);
    //MatDestroy(&submatBT);
  }
//         }
  //PetscFree(colindicesbcpm);
  //PetscFree(colindicesbcpm1);
  printf("B\n");
  if(rank==mpisize-1)proc1=1;
  //ierr = PetscMalloc(proc1*sizeof(Mat *),&submatD);
  //CHKERRQ(ierr);
  printf("D\n");
  MatCreateSubMatrices(submatBB[0],proc1,Bindices,Cindices,MAT_INITIAL_MATRIX,&submatD);
  printf("End partitionning rank %d time %f\n",rank,((double)(clock()-timestr))/CLOCKS_PER_SEC);
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&rowindices[i]);
    CHKERRQ(ierr);
    ierr = ISDestroy(&colindices[i]);
    CHKERRQ(ierr);
    //ierr = ISDestroy(&rowindicesij[i]);
    //CHKERRQ(ierr);
    //ierr = ISDestroy(&colindicesij[i]);
    //CHKERRQ(ierr);
  }
//  for (i=0; i<nmatint; i++) {
//    ierr = ISDestroy(&Cindices[i]);
//    CHKERRQ(ierr);
//    ierr = ISDestroy(&Bindices[i]);
//    CHKERRQ(ierr);
  ierr = ISDestroy(&Cindicesc);
  CHKERRQ(ierr);
  ierr = ISDestroy(&Bindicesc);
  CHKERRQ(ierr);
//  }
  PetscFree(colindices);
  ierr = ISDestroy(&BBindices[0]);
  PetscFree(rowindices);
  //PetscFree(rowindicesij);
  //PetscFree(colindicesij);

  for (i=0; i<nmatint; i++) {
    ierr = ISDestroy(&rowindicesbc[i]);
    CHKERRQ(ierr);
    //ierr = ISDestroy(&colindicesbc[i]);
    //CHKERRQ(ierr);
  }
  PetscFree(rowindicesbc);
  //PetscFree(colindicesbc);

  PetscFree(Cindices);
  PetscFree(Bindices);
  PetscFree(BBindices);
  PetscFree(indexBB);
  PetscFree(indexBBi);
  ierr = MatDestroy(&submatBB[0]);
  PetscFree(submatBB);
  CHKERRQ(ierr);
  printf("OK1 rank %d sumrowcolin %d\n",rank,sumrowcolin);
//  ha_cgetype *xi1 = (ha_cgetype*)calloc(sumrowcolin,sizeof(ha_cgetype));
  long int *bivinzrow=NULL;//(long int *) calloc (1,sizeof(long int));
  PetscInt *bivinzcol=NULL;//(PetscInt *) calloc (1,sizeof(PetscInt));
  PetscReal *ccolnorms=NULL;//(PetscReal *) calloc (1,sizeof(PetscReal));
  //ccolnorms=realloc(ccolnorms,1*sizeof(PetscReal));
  int bivirowsize=1,bivicolsize=1,bbrowij,ddrowi;
  long int vecbivisize=0,li,lj;//,halfvec;
  ha_cgetype *xi1point;
  long int *bivinzrow1=NULL;//(long int *) calloc (1,sizeof(long int));
  long int *bivinzcol1=NULL;//(long int *) calloc (1,sizeof(long int));
  //ha_cgetype *vecbivi= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
  uvadd xi1indx=0;
  long int nz0=0,nz1,nz3;//,nz2,halfj2;
  int fd1,fd2,fd3,frrsl1,frrsl2,frrsl3;
  omp_set_num_threads(smallthreads);
  #pragma omp parallel private(jthrd,nthrd,j4,j1name,filename,frd,fd1,nrow,ncol,i,j,presolfile,nz,nz1,cntl6in,fwrt,ldsize,j3) shared(insize)
  {
  int windx=0,bindx,eindx;
  long int insizeda0=0,insizeda1=0,insizeda2=0;
  jthrd=omp_get_thread_num();
  nthrd=omp_get_max_threads( );
  nthrd=90/nthrd;
  bindx=nthrd*jthrd;
  eindx=bindx+nthrd;
  windx=bindx;
  int *insizeda=(int *) calloc (5+nreg*insizes,sizeof(int));
  #pragma omp for  schedule (static)
  for(j3=0; j3<nmatint; j3++) {
    if(j3<10)strcpy(j1name,"000");
    if(j3<100&&j3>9)strcpy(j1name,"00");
    if(j3<1000&&j3>99)strcpy(j1name,"0");
    if(j3>=1000)j1name[0]='\0';
    sprintf(filename, "%d",j3);
    strcat(j1name,filename);
    strcpy(filename,"_rank");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    //printf("filename %s\n",filename);
    #pragma omp critical
    {
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Cannot open file. Please run presol by setting -presol 1!\n");
    }
    frd=fread(insizeda, sizeof(int), 5+nreg*insizes, presolfile);
    fclose(presolfile);
    if( frd== 0)printf("File read error. No presol!\n");
    }
    if(insizeda0<insizeda[0])insizeda0=insizeda[0];
    if(insizeda1<insizeda[1])insizeda1=insizeda[1];
    if(insizeda2<insizeda[2])insizeda2=insizeda[2];
  }
  ldsize=ceil((laDi/100.0)*insizeda2);
  ldsize+=10;
  insizeda2=ldsize;
  int *irn1=(int *) malloc(ldsize*sizeof(int));
  if(irn1==NULL)printf("Error!!!! No memory left for irn1!\n");
  int *jcn1=(int *) malloc(ldsize*sizeof(int));
  if(jcn1==NULL)printf("Error!!!! No memory left for jcn1!\n");
  ha_cgetype *vecbivi0= (ha_cgetype *) malloc (ldsize*sizeof(ha_cgetype));
  if(vecbivi0==NULL)printf("Error!!!! No memory left for vecbivi0!\n");
  int *irn1a=(int *) malloc (insizeda0*sizeof(int));
  if(irn1a==NULL)printf("Error!!!! No memory left for irn1a!\n");
  int *jcn1a=(int *) malloc (insizeda1*sizeof(int));
  if(jcn1a==NULL)printf("Error!!!! No memory left for jcn1a!\n");
  ha_cgetype *fw=(ha_cgetype *) malloc (insizeda0*sizeof(ha_cgetype));
  if(fw==NULL)printf("Error!!!! No memory left for fw!\n");
  int *fiw=(int *) malloc ((6*insizeda0+3*insizeda1)*sizeof(int));
  if(fiw==NULL)printf("Error!!!! No memory left for fiw!\n");
  int *fkeep=(int *) malloc ((insizeda0+5*insizeda1+4*insizeda1/1+7)*sizeof(int));//ICNTL(6)=1
  if(fkeep==NULL)printf("Error!!!! No memory left for fkeep!\n");
  #pragma omp for  schedule (static)
  for(j3=0; j3<nmatint; j3++) {
    //ha_cgetype *vecbivi0= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
    //xi1indx=0;
    j4=nreg+j3*(nreg+1);
    //int *irn1=(int *) calloc (1,sizeof(int));
    //int *jcn1=(int *) calloc (1,sizeof(int));
    //read from pre
    //int *insizeda=(int *) calloc (5+nreg*insizes,sizeof(int));
    #pragma omp critical
    {
    if(j3<10)strcpy(j1name,"000");
    if(j3<100&&j3>9)strcpy(j1name,"00");
    if(j3<1000&&j3>99)strcpy(j1name,"0");
    if(j3>=1000)j1name[0]='\0';
    sprintf(filename, "%d",j3);
    strcat(j1name,filename);
    strcpy(filename,"_rank");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    //printf("filename %s\n",filename);
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Cannot open file. Please run presol by setting -presol 1!\n");
    }
    frd=fread(insizeda, sizeof(int), 5+nreg*insizes, presolfile);
    fclose(presolfile);
    remove(filename);
    if( frd== 0)printf("File read error. No presol!\n");
    }
    nrow=insizeda[3];
    ncol=insizeda[3];
    #pragma omp critical
    {
    for(i=0; i<nreg; i++) {
      for(j=0; j<insizes; j++) {
        //printf("in1 % in2 %d\n",insize[(i+j3*(nreg+1))*insizes+12],insizeda[i+5]);
        insize[(i+j3*(nreg+1))*insizes+j]=insizeda[i*insizes+5+j];
      }
    }
    //insize[j4*insizes+12]=insizeda[5];
    //printf("in %d %d %d %d\n",insizeda[0],insizeda[1],insizeda[2],insizeda[3]);
    //irn1=realloc(irn1,insizeda[2]*sizeof(int));
    //jcn1=realloc(jcn1,insizeda[2]*sizeof(int));
    //vecbivi0=realloc(vecbivi0,insizeda[2]*sizeof(ha_cgetype));
    //int *irn1a=(int *) calloc (insizeda[0],sizeof(int));
    //int *jcn1a=(int *) calloc (insizeda[1],sizeof(int));
//     #pragma omp critical
//     {
    strcpy(filename,"_row");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Cannot open file.\n");
    }
    frd=fread(irn1a, sizeof(int), insizeda[0], presolfile);
    if(frd == 0)printf("File read error.");
    fclose(presolfile);
    remove(filename);
    strcpy(filename,"_col");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Cannot open file.\n");
    }
    frd=fread(jcn1a, sizeof(int), insizeda[1], presolfile);
    if(frd== 0) printf("File read error.");
    fclose(presolfile);
    remove(filename);
    if(isLinux==0) {
      strcpy(filename,"_bivi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Cannot open file.\n");
      }
      frd=fread(vecbivi0, sizeof(ha_cgetype), insizeda[2], presolfile);
      if(frd== 0) printf("File read error %ld %d.\n",frd,insizeda[2]);
      fclose(presolfile);
      remove(filename);
      strcpy(filename,"_rbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Cannot open file.\n");
      }
      frd=fread(irn1, sizeof(int), insizeda[2], presolfile);
      if(frd== 0) printf("File read error.\n");
      fclose(presolfile);
      remove(filename);
      strcpy(filename,"_cbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Cannot open file.\n");
      }
      frd=fread(jcn1, sizeof(int), insizeda[2], presolfile);
      if(frd== 0) printf("File read error.");
      fclose(presolfile);
      remove(filename);
    }
    else {
      strcpy(filename,"_bivi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      fd1 = open(filename, O_RDONLY);
      frd = read(fd1,vecbivi0, sizeof(ha_cgetype)*insizeda[2]);
      close(fd1);
      remove(filename);
      strcpy(filename,"_rbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      fd1 = open(filename, O_RDONLY);
      frd = read(fd1,irn1, sizeof(int)*insizeda[2]);
      close(fd1);
      remove(filename);
      strcpy(filename,"_cbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      fd1 = open(filename, O_RDONLY);
      frd = read(fd1,jcn1, sizeof(int)*insizeda[2]);
      close(fd1);
      remove(filename);
    }
    }
   
    if(insizeda[0]>insizeda[3]||insizeda[1]>insizeda[3]) {
      patio_mat_(insizeda,irn1,jcn1,vecbivi0,irn1a,jcn1a);
    }
    //free(irn1a);
    //free(jcn1a);
    nz1=insizeda[2];
    //free(insizeda);
    //end read from pre
    //printf("Arank %d j4 %d thrd %d\n",rank,j4,jthrd);
    nz=nz1;
    //obiviindx1=realloc(obiviindx1,sizeof(long int));
    //obiviindx0=realloc(obiviindx0,sizeof(long int));
    ldsize=ceil((laDi/100.0)*nz1);
    //irn1=realloc(irn1,ldsize*sizeof(int));
    //jcn1=realloc(jcn1,ldsize*sizeof(int));
    //vecbivi0=realloc(vecbivi0,ldsize*sizeof(ha_cgetype));
    //memset(vecbivi0+nz1,0,(ldsize-nz1)*sizeof(ha_cgetype));
    //memset(irn1+nz1,0,(laDi-1)*nz1*sizeof(int));
    #pragma omp critical
    {
    insize[j4*insizes]=nrow;
    insize[j4*insizes+1]=ncol;
    insize[j4*insizes+2]=nz1;
    insize[j4*insizes+10]=rank;
    insize[j4*insizes+11]=j4;
    insize[j4*insizes+9]=laDi;
    insize[j4*insizes+16]=ldsize;
    insize[j4*insizes+15]=7+windx;//(j3+jthrd)%90+7;
    }
    windx++;
    if(windx==eindx)windx=bindx;
    //for(i=0;i<nz;i++)printf("aii %f irn1 %d jcn %d\n",vecbivi0[i],irn1[i],jcn[i]);
    //printf("ok!!!!!!!! vecbiuisize %d nrow %d\n",ncol,nrow);
    //printf("ok!!!!!!!! %d %d\n",presol,nrow);
      //prep48_alu_(insize+j4*insizes,irn1,jcn1,vecbivi0);
      prep48_alu1_(insize+j4*insizes,irn1,jcn1,vecbivi0,fw,fiw,fkeep);
    //insize[j4*insizes+15]=0;
    printf("rank %d ok0 nz %ld nzmax %ld ldsize %ld!!!!!!!!\n",rank,nz1,insizeda2,ldsize);
    //vecbivi0=realloc(vecbivi0,1*sizeof(ha_cgetype));
    //free(vecbivi0);
    //free(irn1);
    //free(jcn1);
  }
  printf("rank %d ok1!!!!!!!!\n",rank);
  free(vecbivi0);
  printf("rank %d ok2!!!!!!!!\n",rank);
  free(irn1);
  printf("rank %d ok3!!!!!!!!\n",rank);
  free(jcn1);
  printf("rank %d ok4!!!!!!!!\n",rank);
  free(irn1a);
  printf("rank %d ok5!!!!!!!!\n",rank);
  free(jcn1a);
  printf("rank %d ok6!!!!!!!!\n",rank);
  free(fw);
  printf("rank %d ok7!!!!!!!!\n",rank);
  free(fiw);
  printf("rank %d ok8!!!!!!!!\n",rank);
  free(fkeep);
  printf("rank %d ok9!!!!!!!!\n",rank);
  free(insizeda);
  printf("rank %d ok10!!!!!!!!\n",rank);
  }
  omp_set_num_threads(mymaxnumthrd);
  if(presol)
    return 0;
//  ierr = PetscFree(submatAij);
  printf("Rank %d End of preparation!\n",rank);
  MPI_Barrier(PETSC_COMM_WORLD);

  //bbrowij=submatB[0]->rmap->n;
  //printf("j3 %d bbrowij %d\n",j3,bbrowij);
  ccolnorms=realloc(ccolnorms,BBrow*sizeof(PetscReal));
  bivinzrow=realloc(bivinzrow,BBrow*sizeof(long int));
  memset(bivinzrow,0,BBrow*sizeof(long int));
  bivinzcol=realloc(bivinzcol,BBrow*sizeof(PetscInt));
  memset(bivinzcol,0,BBrow*sizeof(PetscInt));
  printf("rank %d OKKK BBrow %d\n",rank,BBrow);
  for(j1=0; j1<nmatint; j1++) {
    //if(j1<nmatint) {
    //printf("rank %d OKKK1 j1 %d\n",rank,j1);
    MatGetColumnNorms(submatC[j1],NORM_1,ccolnorms);
    //printf("OKKK1 j1 %d ccol %d\n",j1,submatC[j1]->cmap->n);
    //printf("j11111 %d\n",j1);
    for(i=0; i<BBrow; i++)if(ccolnorms[i]>0)bivinzcol[i]++;
    //printf("OKKK1 j1 %d\n",j1);
    Mat_SeqAIJ         *a1b=(Mat_SeqAIJ*)submatB1[j1][0]->data;//*aa=subA->data;
    ai= a1b->i;
    nz=a1b->nz;
    nrow=submatB1[j1][0]->rmap->n;
    for(i=0; i<nrow-1; i++)if(ai[i]!=ai[i+1])bivinzrow[i]++;
    if(ai[i]<nz)bivinzrow[i]++;

    Mat_SeqAIJ         *a2b=(Mat_SeqAIJ*)submatB2[j1][0]->data;//*aa=subA->data;
    ai= a2b->i;
    nz=a2b->nz;
    nrow=submatB2[j1][0]->rmap->n;
    for(j=0; j<nrow-1; j++)if(ai[j]!=ai[j+1])bivinzrow[j]++;
    if(ai[j]<nz)bivinzrow[j]++;
    nz=nz+a1b->nz;
    //nrow=nrow+submatB1[j1][0]->rmap->n;
    //printf("j11111 %d i %d\n",j,i);
    //}
  }
  //printf("OKKK\n");
  free(ccolnorms);
  bivirowsize=0;
  bivicolsize=0;
  //printf("BBrow %d vec-sum %ld\n",BBrow,VecSize-sumrowcolin);
  for(i=0; i<VecSize-sumrowcolin; i++) {
    if(bivinzrow[i]>0) {
      bivinzrow[i]=bivirowsize;
      bivirowsize++;
    }
    else bivinzrow[i]=-1;
    if(bivinzcol[i]>0) {
      bivinzcol[i]=bivicolsize;
      bivicolsize++;
    }
    else bivinzcol[i]=-1;
  }
  //bivirowsize--;
  //bivicolsize--;
  bivinzrow1=realloc (bivinzrow1,bivirowsize*sizeof(long int));//(int *) calloc (bivirowsize,sizeof(int));
  memset(bivinzrow1,0,bivirowsize*sizeof(long int));
  bivinzcol1=realloc (bivinzcol1,bivicolsize*sizeof(long int));//(int *) calloc (bivicolsize,sizeof(int));
  memset(bivinzcol1,0,bivicolsize*sizeof(long int));
  j=0;
  j2=0;
  long int lBBrow=BBrow;
  for(li=0; li<VecSize-sumrowcolin; li++) {
    if(bivinzrow[li]>=0) {
      bivinzrow1[j]=li*lBBrow;
      //printf("j %d bivinzrow1[j] %ld\n",j,bivinzrow1[j]);
      j++;
      bivinzrow[li]=bivinzrow[li]*bivicolsize;//faster MSOL
    }
    if(bivinzcol[li]>=0) {
      bivinzcol1[j2]=li;
      j2++;
    }
  }
  //printf("j %d rows %d j2 %d cols %d\n",j,bivirowsize,j2,bivicolsize);
  vecbivisize=((long int)bivirowsize)*((long int)bivicolsize);
  if(vecbivisize<0)printf("Overflow detected, the interface problem is too big!!!!!!\n");
  ha_cgetype *vecbivi= (ha_cgetype *) calloc (vecbivisize,sizeof(ha_cgetype));
  //vecbivi=realloc(vecbivi,vecbivisize*sizeof(ha_cgetype));
  //memset(vecbivi,0,vecbivisize*sizeof(ha_cgetype));
  printf("vecbivisize %ld rank %d\n",vecbivisize,rank);
  //uvadd xi1indx=0;
  xi1indx=0;
  //j4=0;
  //int *insizeDD=(int *) calloc (insizes*nmatint,sizeof(int));
  //printf("OKt here11 nrowc %d ncolc %d!\n",nrowc,ncolc);
  int maxcolc1=0,maxcolc=0,maxrowc=0,maxrowc1=0;
  long int longsize;
  int *vecnrowc= (int*)calloc(nmatint+1,sizeof(int));
  for(j1=0; j1<nmatint; j1++) {
    MatGetSize(submatC[j1],&maxrowc1,&maxcolc1);
    vecnrowc[j1+1]=vecnrowc[j1]+maxrowc1;
    if(maxcolc<maxcolc1)maxcolc=maxcolc1;
    if(maxcolc<maxrowc1)maxcolc=maxrowc1;
    if(maxrowc<maxrowc1)maxrowc=maxrowc1;
  }
  int maxrowcij;
  ha_cgetype *xi1 = (ha_cgetype*)calloc(sumrowcolin,sizeof(ha_cgetype));
  #pragma omp parallel private(jthrd,timestr,aic,ajc,valsc,nrowc,ncolc,a1i,a1j,val1s,nz,a2i,a2j,val2s,nrowb,ncolb,i,j,j2,xi1point,xi1indx,maxrowcij,la1,fp1,fp2,fp3,freadresult,frrsl1,frrsl2,frrsl3,longsize,nzc) shared(submatC,submatB1,submatB2,xi1,submatCij,submatBij,insize,yi1,vecbivi,vecbiui)
  {
  int icntl[20],info[20];
  ha_cgetype cntl[10],rinfo[10],error1[3];
  ha_cgetype *w= (ha_cgetype*)calloc(4*maxcolc,sizeof(ha_cgetype));
  int *iw= (int*)calloc(maxcolc,sizeof(int));
  int **irnereg= (int**)calloc(nreg+1,sizeof(int*));
  int **keepreg= (int**)calloc(nreg+1,sizeof(int*));
  ha_cgetype **valereg = (ha_cgetype**)calloc(nreg+1,sizeof(ha_cgetype*));
  FILE** nfp1= (FILE**)calloc(nreg+1,sizeof(FILE*));
  FILE** nfp2= (FILE**)calloc(nreg+1,sizeof(FILE*));
  FILE** nfp3= (FILE**)calloc(nreg+1,sizeof(FILE*));
  int* nfd1= (int*)calloc(nreg+1,sizeof(int));
  int* nfd2= (int*)calloc(nreg+1,sizeof(int));
  int* nfd3= (int*)calloc(nreg+1,sizeof(int));
  Mat submatCT;
//   for(i=0; i<nreg+1; i++) {
//     irnereg[i]= (int*)calloc(1,sizeof(int));
//     keepreg[i]= (int*)calloc(1,sizeof(int));
//     valereg[i] = (ha_cgetype*)calloc(1,sizeof(ha_cgetype));
//   }
  jthrd=omp_get_thread_num();
  #pragma omp for
  for(j1=0; j1<nmatint; j1++) {
    xi1indx=vecnrowc[j1];
    timestr=clock();//time(&timestr);//=clock();
    MatTranspose(submatC[j1],MAT_INITIAL_MATRIX,&submatCT);
    Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatCT->data;//*aa=subA->data;
    //insizeDD[j1*insizes+13]=bivirowsize;
    //insizeDD[j1*insizes+14]=bivicolsize;
    aic= ac->i;
    ajc= ac->j;
    valsc=ac->a;
    nzc=ac->nz;
    nrowc=submatCT->rmap->n;
    ncolc=submatCT->cmap->n;
    //insizeDD[j1*insizes+3]=nrowc;
    //insizeDD[j1*insizes+4]=ncolc;
    //insizeDD[j1*insizes+5]=nzc;
    //insizeDD[j1*insizes+10]=rank;
    //insizeDD[j1*insizes+11]=j1;
    //MatGetSize(submatB1[j1][0],&nrowb,&ncolb);
    Mat_SeqAIJ         *a1b=(Mat_SeqAIJ*)submatB1[j1][0]->data;//*aa=subA->data;
    a1i= a1b->i;//irnbs+j1irnbs;
    a1j= a1b->j;//jcnbs+j1nz;
    val1s=a1b->a;//valbs+j1nz;
    nz=a1b->nz;//nzbs[j1];
    //MatGetSize(submatB2[j1][0],&nrowb,&ncolb);
    Mat_SeqAIJ         *a2b=(Mat_SeqAIJ*)submatB2[j1][0]->data;//*aa=subA->data;
    a2i= a2b->i;//irnbs+j1irnbs;
    a2j= a2b->j;//jcnbs+j1nz;
    val2s=a2b->a;//valbs+j1nz;
    nz=nz+a2b->nz;//nzbs[j1];
    nrowb=submatB1[j1][0]->rmap->n;
    ncolb=submatB1[j1][0]->cmap->n;
    PetscInt *ai12= (PetscInt*)calloc(nrowb+1,sizeof(PetscInt));
    PetscInt *aj12= (PetscInt*)calloc(nz,sizeof(PetscInt));
    PetscScalar *vals12= (PetscScalar*)calloc(nz,sizeof(PetscScalar));
    for(i=0; i<nrowb-1; i++) {
      ai12[i]=a1i[i]+a2i[i];
      memcpy(&aj12[ai12[i]],&a1j[a1i[i]],(a1i[i+1]-a1i[i])*sizeof(PetscInt));
      memcpy(&vals12[ai12[i]],&val1s[a1i[i]],(a1i[i+1]-a1i[i])*sizeof(PetscScalar));
      j2=ai12[i]+a1i[i+1]-a1i[i];
      memcpy(&vals12[j2],&val2s[a2i[i]],(a2i[i+1]-a2i[i])*sizeof(PetscScalar));
      for(j=a2i[i]; j<a2i[i+1]; j++)aj12[j2+j-a2i[i]]=a2j[j]+ncolb;
    }
    ai12[i]=a1i[i]+a2i[i];
    ai12[i+1]=nz;
    nz=a1b->nz;
    memcpy(&aj12[ai12[i]],&a1j[a1i[i]],(nz-a1i[i])*sizeof(PetscInt));
    memcpy(&vals12[ai12[i]],&val1s[a1i[i]],(nz-a1i[i])*sizeof(PetscScalar));
    j2=ai12[i]+nz-a1i[i];
    nz=a2b->nz;
    memcpy(&vals12[j2],&val2s[a2i[i]],(nz-a2i[i])*sizeof(PetscScalar));
    for(j=a2i[i]; j<nz; j++)aj12[j2+j-a2i[i]]=a2j[j]+ncolb;
    nz=nz+a1b->nz;
    ncolb=ncolb+submatB2[j1][0]->cmap->n;
    MatDestroy(&submatB1[j1][0]);//submatBT);
    free(submatB1[j1]);
    MatDestroy(&submatB2[j1][0]);//submatBT);
    free(submatB2[j1]);
    //insizeDD[j1*insizes+6]=nrowb;//ncolb;
    //insizeDD[j1*insizes+7]=ncolb;//nrowb;
    //insizeDD[j1*insizes+8]=nz;
    //for(i=0;i<nz;i++)++aj[i];
    xi1point=xi1+xi1indx;
    printf("OK here j1 %d yi1 %ld xi1indx %ld a2j %d!\n",j1,sizeof(yi1[j1])/sizeof(PetscInt),xi1indx,ncolb);
    maxrowcij=0;
    #pragma omp critical
    {
    for(i=0; i<nreg+1; i++) {
      j2=j1*(nreg+1)+i;
      insize[j2*insizes+15]=j2%90+7;
      if(i!=nreg&&(submatCij[j2][0]->rmap->n)>maxrowcij)maxrowcij=submatCij[j2][0]->rmap->n;
      la1=ceil((insize[j2*insizes+9]/100.0)*insize[j2*insizes+2]);
      if(isLinux==0) {
        nfp1[i] = fopen(fn01[j2], "rb");
        //printf("fn1 %s\n",fn01[j2]);
        if (nfp1[i]==NULL)printf("File opening error\n");
        irnereg[i] = realloc(irnereg[i],la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
        freadresult=fread(irnereg[i],sizeof(int),la1,nfp1[i]);
        fclose(nfp1[i]);
        //printf("fn2 %s\n",fn01[j2]);
        nfp2[i] = fopen(fn02[j2], "rb");
        if (nfp2[i]==NULL)printf("File opening error\n");
        keepreg[i] = realloc(keepreg[i],insize[j2*insizes+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
        freadresult=fread(keepreg[i],sizeof(int),insize[j2*insizes+12],nfp2[i]);
        fclose(nfp2[i]);
        //printf("fn3 %s\n",fn01[j2]);
        nfp3[i] = fopen(fn03[j2], "rb");
        if (nfp3[i]==NULL)printf("File opening error\n");
        valereg[i] = realloc(valereg[i],la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
        freadresult=fread(valereg[i],sizeof(ha_cgetype),la1,nfp3[i]);
        fclose(nfp3[i]);
        //printf("fn4 %s\n",fn01[j2]);
      }
      else {
        nfd1[i] = open(fn01[j2],O_RDONLY);
        nfd2[i] = open(fn02[j2],O_RDONLY);
        nfd3[i] = open(fn03[j2],O_RDONLY);
        //la1=insize[j2*insizes+9]*insize[j2*insizes+2];
        la1=ceil((insize[j2*insizes+9]/100.0)*insize[j2*insizes+2]);
        irnereg[i] = realloc(irnereg[i],la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
        keepreg[i] = realloc(keepreg[i],insize[j2*insizes+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
        valereg[i] = realloc(valereg[i],la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
        frrsl1=read(nfd1[i],irnereg[i],sizeof(int)*la1);
        frrsl2=read(nfd2[i],keepreg[i],insize[j2*insizes+12]*sizeof(int));
        frrsl3=read(nfd3[i],valereg[i],la1*sizeof(ha_cgetype));
        close(nfd1[i]);
        close(nfd2[i]);
        close(nfd3[i]);
      }
    }
    }
    printf("OK here11a nrowc %d ncolc %d maxrow %d!\n",nrowc,ncolc,maxrowcij);
    ha_cgetype *b02 = (ha_cgetype*)calloc(maxrowcij,sizeof(ha_cgetype));

    NDBBD_sol_nread1(rank,j1*(nreg+1),nreg,insize,insizes,submatCij,submatBij,yi1[j1],xi1point,irnereg,keepreg,valereg,cntl,rinfo,error1,icntl,info,w,iw,b02);
    //printf("Submatrix %d rank %d one mat calculation time %f\n",j1,rank,difftime(timeend,timestr));
    //if(j1==0)for(i=1;i<1000;i++)printf("j1 %d i %d xi %lf\n",j1,i,xi1point[i]);
    //printf("OK here11b nrowc %d ncolc %d rank %d!\n",nrowc,ncolc,rank);
    ha_cgetype *bccol= (ha_cgetype*)calloc(ncolc,sizeof(ha_cgetype));
    longsize=ncolc*sizeof(ha_cgetype);
    for(i=0; i<nrowc-1; i++) {
      if(aic[i]<aic[i+1]) {
        memset(yi1[j1],0,longsize);
        for(j=aic[i]; j<aic[i+1]; j++) {
          yi1[j1][ajc[j]]=valsc[j];
        }
        //printf("OK here11 i %d nrowc %d!\n",i,nrowc);
        NDBBD_sol_nread1(rank,j1*(nreg+1),nreg,insize,insizes,submatCij,submatBij,yi1[j1],bccol,irnereg,keepreg,valereg,cntl,rinfo,error1,icntl,info,w,iw,b02);
        //printf("OK here112 i %d!\n",i);
        spar_vbiviadd_(bccol,bivinzcol+i,bivinzrow,&vecbivisize,&nrowb,&ncolc,&nz,ai12,aj12,vals12,vecbivi);
//        for(j2=0; j2<nrowb-1; j2++) {
//          for(j=ai[j2]; j<ai[j2+1]; j++) {
//            vecbivi[bivinzcol[i]+bivinzrow[j2]]-=vals[j]*bccol[aj[j]];//xi[i];
//            //printf("rank %d j1 %d j %d bccol %f",rank,j1,j,bccol[aj[j]]);
//          }
//        }
//        //printf("OK here112 i %d!\n",i);
//        for(j=ai[j2]; j<nz; j++) {
//          vecbivi[bivinzcol[i]+bivinzrow[j2]]-=vals[j]*bccol[aj[j]];//xi[i];
//          //printf("rank %d j1 %d j %d bccol %f",rank,j1,j,bccol[aj[j]]);
//        }
      }
    }
    //printf("OK here11c nrowc %d ncolc %d rank %d!\n",nrowc,ncolc,rank);
    if(aic[i]<nzc) {
      memset(yi1[j1],0,longsize);
      for(j=aic[i]; j<nzc; j++) {
        yi1[j1][ajc[j]]=valsc[j];
      }
      NDBBD_sol_nread1(rank,j1*(nreg+1),nreg,insize,insizes,submatCij,submatBij,yi1[j1],bccol,irnereg,keepreg,valereg,cntl,rinfo,error1,icntl,info,w,iw,b02);
      spar_vbiviadd_(bccol,bivinzcol+i,bivinzrow,&vecbivisize,&nrowb,&ncolc,&nz,ai12,aj12,vals12,vecbivi);
//      for(j2=0; j2<nrowb-1; j2++) {
//        for(j=ai[j2]; j<ai[j2+1]; j++) {
//          vecbivi[bivinzcol[i]+bivinzrow[j2]]-=vals[j]*bccol[aj[j]];//xi[i];
//        }
//      }
//      for(j=ai[j2]; j<nz; j++) {
//        vecbivi[bivinzcol[i]+bivinzrow[j2]]-=vals[j]*bccol[aj[j]];//xi[i];
//      }
    }
    //printf("OK here11d nrowc %d ncolc %d rank %d!\n",nrowc,ncolc,rank);
    free(bccol);
    free(b02);
    //spec48_msol_(insize+j1*insizes,irn,jcn,values,yi1[j1],xi1point,aic,ajc,valsc,ai,aj,vals,vecbivi,bivinzrow,bivinzcol);
    MatDestroy(&submatCT);
    free(yi1[j1]);
    //Multiply Bi by ui:
    spar_mulmin_(xi1point,&nrowb,&nz,ai12,aj12,vals12,vecbiui);
//    for(i=0; i<nrowb-1; i++) {
//      for(j=ai[i]; j<ai[i+1]; j++) {
//        vecbiui[i]-=vals[j]*xi1point[aj[j]];//xi[i];
//      }
//    }
//    for(j=ai[i]; j<nz; j++) {
//      vecbiui[i]-=vals[j]*xi1point[aj[j]];//xi[i];
//      //printf("i %d ui %lf\n",aj[j],xi1point[aj[j]]);
//    }
    //if(j1==10)for(j=0; j<ncolc; j++)printf("rank %d j1 %d x1 %f\n",rank,j1,xi1point[j]);
    free(ai12);
    free(aj12);
    free(vals12);
    //xi1indx+=ncolc;
//    for(i=0; i<nreg; i++) {
//      //MatDestroy(&submatBij[j1*(nreg+1)+i][0]);
//      MatDestroy(&submatCij[j1*(nreg+1)+i][0]);
//    }
    timeend=clock();//(&timeend);
    printf("Submatrix %d rank %d calculation time %f\n",j1,rank,difftime(timeend,timestr)/CLOCKS_PER_SEC);
    //j4+=ncolc;
  }
  free(w);
  free(iw);
  for(i=0; i<nreg+1; i++) {
    free(irnereg[i]);//= (int*)calloc(1,sizeof(int));
    free(keepreg[i]);//= (int*)calloc(1,sizeof(int));
    free(valereg[i]);// = (ha_cgetype*)calloc(1,sizeof(ha_cgetype));
  }
  free(irnereg);//= (int**)calloc(nreg+1,sizeof(int*));
  free(keepreg);//= (int**)calloc(nreg+1,sizeof(int*));
  free(valereg);// = (ha_cgetype**)calloc(nreg+1,sizeof(ha_cgetype*));
  free(nfp1);
  free(nfp2);
  free(nfp3);
  free(nfd1);
  free(nfd2);
  free(nfd3);
  }
  printf("Completed calculation of partitioned matrices! Rank %d\n",rank);
  free(bivinzcol);
  free(bivinzrow);
  free(submatB1);
  free(submatB2);
  //ierr = PetscFree(submatBij);
//  for(i=0; i<nmatin; i++)free(submatCij[i]);
//  free(submatCij);
  CHKERRQ(ierr);
  free(yi1);
  //MPI_Barrier(PETSC_COMM_WORLD);
  printf("Completed calculation of partitioned matrices! Rank %d\n",rank);
  timestr=clock();

  #pragma omp parallel private(li) reduction(+:lj)
  {
  lj=0;
  #pragma omp for
  //for(li=vecbivisize-1; li!=-1; --li){
  for(li=0; li<vecbivisize; ++li){
    if(vecbivi[li]!=0)lj++;
  }
  }
  obiviindx1=realloc(obiviindx1,(lj+1)*sizeof(long int));
  obiviindx1[0]=-1;
  biviindx1=obiviindx1;
  biviindx1++;
  //int nz0,nz1,nz2,nz3=j,halfj2;
  long int lj2=0;
  printf("vbivi %ld colsize %d\n",vecbivisize,bivicolsize);
  for(li=0; li<vecbivisize; ++li) {
    if(vecbivi[li]!=0) {
      vecbivi[lj2]=vecbivi[li];
      biviindx1[lj2]=bivinzcol1[li%bivicolsize]+bivinzrow1[li/bivicolsize];//[i-j1*bivicolsize]+bivinzrow1[j1];//*nrowb;
      //if(rank==1)printf("j2 %d bivi %ld c %d r %ld\n",j2,biviindx1[j2],bivinzcol1[i%bivicolsize],bivinzrow1[i/bivicolsize]);
      //if(biviindx1[j2]>=vecbivisize)printf("b %d v %d\n",biviindx1[j2],vecbivisize);
      ++lj2;
    }
  }
  //for(i=0; i<bivirowsize; ++i)printf("OK here11! i %d j %d j2 %d\n",i,bivinzrow1[i],bivinzcol1[i]);
  free(bivinzrow1);
  free(bivinzcol1);
  printf("OK here11! j %ld j2 %ld bivi %ld\n",lj,lj2,biviindx1[lj2-1]);
  vecbivi=realloc(vecbivi,lj2*sizeof(ha_cgetype));
  printf("OK here11!\n");
  nz0=lj2;
  int group_size,ha_id,group_size1,group_size11,ha_id1,color;
  MPI_Comm_rank( HA_COMM, &ha_id);
  MPI_Comm_size(HA_COMM,&group_size);
  //printf("rank %d color %d group %d\n",rank,color,group_size);
  if(ha_id==group_size-1)color=1;
  else color=0;
  MPI_Comm_rank( HA1_COMM, &ha_id1);
  MPI_Comm_size(HA1_COMM,&group_size1);
  group_size11=group_size1;
  MPI_Bcast(&group_size11,1, MPI_INT,mpisize-1, PETSC_COMM_WORLD);
  //printf("rank %d color %d group11 %d\n",rank,color,group_size11);
  //printf("rank %d j2 %d!!!!!!!!!!!!\n",rank,j2);
  long int sindx01,sindx02,sindx03;
  
  if(ha_id!=group_size-1) {
    //j2=nz3;
    MPI_Send(&lj2,1, MPI_LONG,group_size-1, 10, HA_COMM);
//    if(SORD==1)MPI_Send(vecbivi,j2, MPI_DOUBLE,group_size-1, 11, HA_COMM);
//    else MPI_Send(vecbivi,j2, MPI_FLOAT,group_size-1, 11, HA_COMM);
//    MPI_Send(biviindx1,j2, MPI_LONG,group_size-1, 12, HA_COMM);
    sindx01=lj2/MAXSSIZE;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
    if(SORD==1)MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_DOUBLE,group_size-1, 11*sindx01+sindx02, HA_COMM);
    else MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_FLOAT,group_size-1, 11*sindx01+sindx02, HA_COMM);
    sindx03+=MAXSSIZE;
    }
    if(SORD==1)MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,group_size-1, 11, HA_COMM);
    else MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,group_size-1, 11, HA_COMM);
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
    MPI_Send(biviindx1+sindx03,MAXSSIZE, MPI_LONG,group_size-1, 12*sindx01+sindx02, HA_COMM);
    sindx03+=MAXSSIZE;
    }
    MPI_Send(biviindx1+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,group_size-1, 12, HA_COMM);
    free(obiviindx1);
    obiviindx1=NULL;
    //obiviindx1=realloc(obiviindx1,sizeof(long int));
    free(vecbivi);
    vecbivi=NULL;
    //vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
  }
  if(ha_id==group_size-1) {
    for(j6=0; j6<group_size-1; j6++) {
      MPI_Recv(&lj2,1, MPI_LONG,  MPI_ANY_SOURCE,10, HA_COMM,&status);
      i=status.MPI_SOURCE;
      //printf("rank1 %d j2 %d!!!!!!!!!!!!\n",rank,j2);
      obiviindx0=realloc(obiviindx0,(lj2+1)*sizeof(long int));
      obiviindx0[0]=-1;
      biviindx0=obiviindx0;
      biviindx0++;
      //printf("rank2 %d j2 %d!!!!!!!!!!!!\n",rank,j2);
      ha_cgetype *vecbivi0= (ha_cgetype *) calloc (lj2,sizeof(ha_cgetype));
//      if(SORD==1)MPI_Recv(vecbivi0,j2, MPI_DOUBLE,i,11, HA_COMM,&status);
//      else MPI_Recv(vecbivi0,j2, MPI_FLOAT,i,11, HA_COMM,&status);
//      MPI_Recv(biviindx0,j2, MPI_LONG,i,12, HA_COMM,&status);
    sindx01=lj2/MAXSSIZE;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      if(SORD==1)MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_DOUBLE,i,11*sindx01+sindx02, HA_COMM,&status);
      else MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_FLOAT,i,11*sindx01+sindx02, HA_COMM,&status);
      sindx03+=MAXSSIZE;
    }
      if(SORD==1)MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,i,11, HA_COMM,&status);
      else MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,i,11, HA_COMM,&status);
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      MPI_Recv(biviindx0+sindx03,MAXSSIZE, MPI_LONG,i,12*sindx01+sindx02, HA_COMM,&status);
      sindx03+=MAXSSIZE;
    }
      MPI_Recv(biviindx0+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,i,12, HA_COMM,&status);
      //printf("rank3 %d j2 %d!!!!!!!!!!!!\n",rank,j2);
      my_spar_compl_(biviindx1,&nz0,biviindx0,&lj2,&nz1);
      //printf("rank4 %d nz0 %d j2 %d nz1 %d!!!!!!!!!!!!\n",rank,nz0,j2,nz1);
      vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
      obiviindx1=realloc(obiviindx1,(nz1+1)*sizeof(long int));
      obiviindx1[0]=-1;
      biviindx1=obiviindx1;
      biviindx1++;
      //printf("rank5 %d j2 %d!!!!!!!!!!!!\n",rank,j2);
      my_spar_add3l_(vecbivi,obiviindx1,&nz0,vecbivi0,obiviindx0,&lj2,&nz1);
      //printf("rank6 %d j2 %d!!!!!!!!!!!!\n",rank,j2);
      free(obiviindx0);
      obiviindx0=NULL;
      //obiviindx0=realloc(obiviindx0,sizeof(long int));
      //free(obiviindx0);
      free(vecbivi0);
      vecbivi0=NULL;
      printf("!!!!!!!!!!!! rank11 %d nz0 %ld j2 %ld nz1 %ld\n",rank,nz0,lj2,nz1);
      nz0=nz1;
    }
  }

  if(group_size11>1&&color==1) {
    nz3=nz0;
    if(ha_id1!=group_size1-1) {
      lj2=nz3;
      MPI_Send(&lj2,1, MPI_LONG,group_size1-1, 10, HA1_COMM);
//      if(SORD==1)MPI_Send(vecbivi,j2, MPI_DOUBLE,group_size1-1, 11, HA1_COMM);
//      else MPI_Send(vecbivi,j2, MPI_FLOAT,group_size1-1, 11, HA1_COMM);
//      MPI_Send(biviindx1,j2, MPI_LONG,group_size1-1, 12, HA1_COMM);
    sindx01=lj2/MAXSSIZE;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      if(SORD==1)MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_DOUBLE,group_size1-1, 11*sindx01+sindx02, HA1_COMM);
      else MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_FLOAT,group_size1-1, 11*sindx01+sindx02, HA1_COMM);
      sindx03+=MAXSSIZE;
    }
      if(SORD==1)MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,group_size1-1, 11, HA1_COMM);
      else MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,group_size1-1, 11, HA1_COMM);
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      MPI_Send(biviindx1+sindx03,MAXSSIZE, MPI_LONG,group_size1-1, 12*sindx01+sindx02, HA1_COMM);
      sindx03+=MAXSSIZE;
    }
      MPI_Send(biviindx1+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,group_size1-1, 12, HA1_COMM);
      free(obiviindx1);
      obiviindx1=NULL;
      //obiviindx1=realloc(obiviindx1,sizeof(long int));
      free(vecbivi);
      vecbivi=NULL;
      //vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
    }
    if(ha_id1==group_size1-1) {
      for(j6=0; j6<group_size1-1; j6++) {
        MPI_Recv(&lj2,1, MPI_LONG,  MPI_ANY_SOURCE,10, HA1_COMM,&status);
        i=status.MPI_SOURCE;
        obiviindx0=realloc(obiviindx0,(lj2+1)*sizeof(long int));
        obiviindx0[0]=-1;
        biviindx0=obiviindx0;
        biviindx0++;
        ha_cgetype *vecbivi0= (ha_cgetype *) calloc (lj2,sizeof(ha_cgetype));
//        if(SORD==1)MPI_Recv(vecbivi0,j2, MPI_DOUBLE,i,11, HA1_COMM,&status);
//        else MPI_Recv(vecbivi0,j2, MPI_FLOAT,i,11, HA1_COMM,&status);
//        MPI_Recv(biviindx0,j2, MPI_LONG,i,12, HA1_COMM,&status);
    sindx01=lj2/MAXSSIZE;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
        if(SORD==1)MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_DOUBLE,i,11*sindx01+sindx02, HA1_COMM,&status);
        else MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_FLOAT,i,11*sindx01+sindx02, HA1_COMM,&status);
        sindx03+=MAXSSIZE;
    }
        if(SORD==1)MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,i,11, HA1_COMM,&status);
        else MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,i,11, HA1_COMM,&status);
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
        MPI_Recv(biviindx0+sindx03,MAXSSIZE, MPI_LONG,i,12*sindx01+sindx02, HA1_COMM,&status);
        sindx03+=MAXSSIZE;
    }
        MPI_Recv(biviindx0+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,i,12, HA1_COMM,&status);
        my_spar_compl_(biviindx1,&nz0,biviindx0,&lj2,&nz1);
        vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
        obiviindx1=realloc(obiviindx1,(nz1+1)*sizeof(long int));
        obiviindx1[0]=-1;
        biviindx1=obiviindx1;
        biviindx1++;
        my_spar_add3l_(vecbivi,obiviindx1,&nz0,vecbivi0,obiviindx0,&lj2,&nz1);
        free(obiviindx0);
        obiviindx0=NULL;
        //obiviindx0=realloc(obiviindx0,sizeof(long int));
        free(vecbivi0);
        vecbivi0=NULL;
        printf("!!!!!!!!!!!! rank22 %d nz0 %ld j2 %ld nz1 %ld\n",rank,nz0,lj2,nz1);
        nz0=nz1;
      }
    }
  }
  forint vecbiuisize=0;
  vecbiuisize=VecSize-sumrowcolin;
  printf("after bivi rank %d\n",rank);
  HaReduceNoComp(vecbiui,vecbiuisize,mpisize,rank,mpisize-1);
  printf("Completed MPI_Reduce Operation! Rank %d time %f\n",rank,((double)(clock()-timestr))/CLOCKS_PER_SEC);
  if(rank!=mpisize-1){
    free(vecbiui);
    vecbiui=NULL;
    //vecbiui=realloc(vecbiui,sizeof(ha_cgetype));//free(vecbiui);
  }
  timestr=clock();
  ha_cgetype *xd;//= (ha_cgetype*)calloc(vecbiuisize,sizeof(ha_cgetype));
  ha_cgetype *x0;//,*vecbivi0;//= (ha_cgetype*)calloc(VecSize,sizeof(ha_cgetype));
  //int j7;
  long int lnz;
  if(rank==mpisize-1) {
    Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatD[0]->data;//*aa=subA->data;
    ai= aa->i;
    aj= aa->j;
    vals=aa->a;
    lnz=aa->nz;
    nrow=submatD[0]->rmap->n;
    ncol=submatD[0]->cmap->n;
    printf("D ncol %d nrow %d nz %ld\n",ncol,nrow,lnz);
    obiviindx0=realloc(obiviindx0,(lnz+1)*sizeof(long int));
    obiviindx0[0]=-1;
    biviindx0=obiviindx0;
    biviindx0++;
    lj2=0;
    lBBrow=ncol;
    for(li=0; li<nrow-1; li++) {
      j7l=li*lBBrow;
      for(j=ai[li]; j<ai[li+1]; j++) {
        biviindx0[lj2]=aj[j]+j7l;
        lj2++;
      }
    }
    j7l=li*lBBrow;
    for(j=ai[li]; j<lnz; j++) {
      biviindx0[lj2]=aj[j]+j7l;
      lj2++;
    }
    my_spar_compl_(biviindx1,&nz0,biviindx0,&lnz,&nz1);
    printf("nz0 %ld nz %ld nz1 %ld\n",nz0,lnz,nz1);
    ldsize=ceil((laD/100.0)*nz1);
    vecbivi=realloc(vecbivi,ldsize*sizeof(ha_cgetype));
    int *irn1=(int *) calloc (nz1,sizeof(int));
    int *jcn=(int *) calloc (ldsize,sizeof(int));
    printf("nz1 %ld obin %ld ond1 %ld\n",nz1,obiviindx1[0],obiviindx0[0]);
    my_spar_add4l_(vecbivi,obiviindx1,irn1,jcn,&nz0,vals,obiviindx0,&lnz,&nz1,&ncol);
    printf("nz1 %ld\n",nz1);
    MatDestroy(&submatD[0]);
    lnz=nz1;
    free(obiviindx1);
    obiviindx1=NULL;
    //obiviindx1=realloc(obiviindx1,sizeof(long int));
    free(obiviindx0);
    obiviindx0=NULL;
    //obiviindx0=realloc(obiviindx0,sizeof(long int));
    irn1=realloc(irn1,ldsize*sizeof(int));
    int *insizeD=(int *) calloc (6,sizeof(int));
    insizeD[0]=nrow;
    insizeD[1]=ncol;
    insizeD[2]=lnz;
    insizeD[3]=laD;
    insizeD[5]=ldsize;
    //printf("ok\n");
    xd=(ha_cgetype *) calloc (vecbiuisize,sizeof(ha_cgetype));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
    //printf("ok!!!!!!!! vecbiuisize %ld nrow %d\n",vecbiuisize,nrow);
    spec48m_ssol2la_(insizeD,irn1,jcn,vecbivi,vecbiui,xd);
    //printf("ok!!!!!!!!\n");
    free(vecbivi);
    vecbivi=NULL;
    //vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
    free(vecbiui);
    vecbiui=NULL;
    //vecbiui=realloc(vecbiui,1*sizeof(ha_cgetype));
    //free(yd);
    free(irn1);
    free(jcn);
    //free(values);
    free(insizeD);
    x0=(ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));//realloc (x0,*sizeof(ha_cgetype));
    for(i=0; i<vecbiuisize; i++) {
      //VecGetValues(xsol,1,&i,&xd[i]);
      x0[indicesC[i]]=xd[i];
      //printf("i %d x1d %lf\n",i,xd[i]);
    }
    //VecView(xsol,0);
    //ierr = VecDestroy(&xsol);
    //CHKERRQ(ierr);
  }
  free(vecbiui);
  obiviindx1=realloc(obiviindx1,sizeof(long int));
  free(obiviindx1);
  obiviindx0=realloc(obiviindx0,sizeof(long int));
  free(obiviindx0);
  free(vecbivi);
  free(indicesC);
  if(rank!=mpisize-1) {
    xd=(ha_cgetype *) calloc (vecbiuisize,sizeof(ha_cgetype));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
    x0=(ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));//realloc (x0,*sizeof(ha_cgetype));
  }
  PetscFree(submatD);
  //PetscFree(submatCT);
  MPI_Barrier(PETSC_COMM_WORLD);
  if(SORD==1)MPI_Bcast(xd, vecbiuisize, MPI_DOUBLE,mpisize-1, PETSC_COMM_WORLD);
  else MPI_Bcast(xd, vecbiuisize, MPI_FLOAT,mpisize-1, PETSC_COMM_WORLD);
  printf("Interface problem solution Rank %d time %f\n",rank,((double)(clock()-timestr))/CLOCKS_PER_SEC);
  timestr=clock();
  xi1indx=0;

  omp_set_num_threads(smallthreads);
  #pragma omp parallel private(jthrd,j1,j2,j,i,timestr,ai,aj,vals,nz,nrow,xi1point,xi1indx,ierr,nindices) shared(submatC,submatCij,submatBij,insize,xi1,colindicesbc1,colindicesbc2,x0)
  {
  jthrd=omp_get_thread_num();
  ha_cgetype *biui0= (ha_cgetype *) calloc (maxrowc,sizeof(ha_cgetype));
  ha_cgetype *be0 = (ha_cgetype*)calloc(maxrowc,sizeof(ha_cgetype));
  #pragma omp for
  for(j1=0; j1<nmatint; j1++) {
    xi1indx=vecnrowc[j1];
    timestr=clock();
    j2=j1*(nreg+1);
    Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatC[j1]->data;//*aa=subA->data;
    //printf("0 rank %d j1 %d vec %lf\n",rank,j1,vecbivi[22800]);
    ai= ac->i;
    aj= ac->j;
    vals=ac->a;
    nz=ac->nz;
    nrow=submatC[j1]->rmap->n;
    //ncol=submatC[j1]->cmap->n;
    //printf("nrow %d\n",nrow);
//    ha_cgetype *biui0= (ha_cgetype *) calloc (nrow,sizeof(ha_cgetype));
//    ha_cgetype *be0 = (ha_cgetype*)calloc(nrow,sizeof(ha_cgetype));
    spar_mulnoadd_(xd,&nrow,&nz,ai,aj,vals,be0);
//    for(i=0; i<nrow-1; i++) {
//      for(j=ai[i]; j<ai[i+1]; j++) {
//        be0[i]+=vals[j]*xd[aj[j]];
//      }
//    }
//    for(j=ai[i]; j<nz; j++) {
//      be0[i]+=vals[j]*xd[aj[j]];
//    }
    MatDestroy(&submatC[j1]);
    //printf("OK here j1 %d yi1 %ld xi1indx %ld!\n",j1,sizeof(yi1[j1])/sizeof(PetscInt),xi1indx);
    ifremove=true;
    NDBBD_sol(rank,j2,nreg,insize,insizes,submatCij,submatBij,be0,biui0,ifremove,fn01,fn02,fn03);
    //printf("OK here11 nrowc %d ncolc %d!\n",nrow,ncol);
//    free(be0);
    //MatDestroy(&submatB[j1]);//submatBT);
    //xi1indx+=nrow;
    for(i=0; i<nreg; i++) {
      MatDestroy(&submatCij[j2+i][0]);
      MatDestroy(&submatBij[j2+i][0]);
    }
    xi1point=xi1+xi1indx;
    ISGetIndices(colindicesbc1[j1],&nindices);
    ISGetSize(colindicesbc1[j1],&j);
    for(i=0; i<j; i++) {
      x0[nindices[i]]+=xi1point[i]-biui0[i];//vals[i]-biui0[i];//vecval;
    }
    ISRestoreIndices(colindicesbc1[j1],&nindices);
    ierr = ISDestroy(&colindicesbc1[j1]);
    ISGetIndices(colindicesbc2[j1],&nindices);
    for(j=i; j<nrow; j++) {
      x0[nindices[j-i]]+=xi1point[j]-biui0[j];//vals[i]-biui0[i];//vecval;
    }
    ISRestoreIndices(colindicesbc2[j1],&nindices);
    ierr = ISDestroy(&colindicesbc2[j1]);

    //CHKERRQ(ierr);
//    free(biui0);
    //xi1indx+=nrow;
    timeend=clock();//time(&timeend);
    printf("Submatrix %d rank %d calculation time %f\n",j1,rank,difftime(timeend,timestr)/CLOCKS_PER_SEC);
    //j4+=ncolc;
  }
  free(be0);
  free(biui0);
  }
  omp_set_num_threads(mymaxnumthrd);
  free(vecnrowc);
  for(i=0; i<nmatin; i++) {
    ierr = PetscFree(submatCij[i]);
    ierr = PetscFree(submatBij[i]);
  }
  ierr = PetscFree(colindicesbc1);
  ierr = PetscFree(colindicesbc2);
  free(submatCij);
  free(submatBij);
  for(j1=0; j1<nmatinplus; j1++) {
    free(fn01[j1]);
    free(fn02[j1]);
    free(fn03[j1]);
  }
  free(fn01);
  free(fn02);
  free(fn03);

//  for(j1=0; j1<nmatint; j1++) {
//    //if(j1<nmatin) {
//      ha_cgetype *biui0= (ha_cgetype *) calloc (ha_ndblocks[j1+begblock[rank]],sizeof(ha_cgetype));
//      fp1 = fopen(fn01[j1], "rb");
//      //printf("fn1 %s\n",fn01[j1]);
//      if (fp1==NULL)printf("File opening error\n");
//      fp2 = fopen(fn02[j1], "rb");
//      if (fp2==NULL)printf("File opening error\n");
//      fp3 = fopen(fn03[j1], "rb");
//      if (fp3==NULL)printf("File opening error\n");
//      la1=insize[j1*insizes+9]*insize[j1*insizes+2];
//      int *irne = (int*)calloc(la1,sizeof(int));
//      int *keep = (int*)calloc(insize[j1*insizes+12],sizeof(int));
//      ha_cgetype *vale = (ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
//      fread(irne,sizeof(int),la1,fp1);
//      fread(keep,sizeof(int),insize[j1*insizes+12],fp2);
//      fread(vale,sizeof(ha_cgetype),la1,fp3);
//      fclose(fp1);
//      fclose(fp2);
//      fclose(fp3);
//      remove(fn01[j1]);
//      remove(fn02[j1]);
//      remove(fn03[j1]);
//      //for(i=0;i<la;i++)printf("vav %f\n",vale[i]);
//      Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatC[j1]->data;//*aa=subA->data;
//      //printf("0 rank %d j1 %d vec %lf\n",rank,j1,vecbivi[22800]);
//      ai= ac->i;
//      aj= ac->j;
//      vals=ac->a;
//      nz=ac->nz;
//      nrow=submatC[j1]->rmap->n;
//      ncol=submatC[j1]->cmap->n;
//      //printf("nrow %d\n",nrow);
//      ha_cgetype *be0 = (ha_cgetype*)calloc(nrow,sizeof(ha_cgetype));
//      for(i=0; i<nrow-1; i++) {
//        for(j=ai[i]; j<ai[i+1]; j++) {
//          be0[i]+=vals[j]*xd[aj[j]];
//        }
//      }
//      for(j=ai[i]; j<nz; j++) {
//        be0[i]+=vals[j]*xd[aj[j]];
//      }
//      MatDestroy(&submatC[j1]);
//      spec48_esol_(insize+j1*insizes,irne,vale,keep,be0,biui0);
//      //printf("here!\n");
//      free(irne);
//      free(keep);
//      free(vale);
//      free(be0);
//      ISGetIndices(colindices[j1],&nindices);
//      xi1point=xi1+xi1indx;
//      for(i=0; i<ha_ndblocks[j1+begblock[rank]]; i++) {
//        x0[nindices[i]]+=xi1point[i]-biui0[i];//vals[i]-biui0[i];//vecval;
//      }
//      ISRestoreIndices(colindices[j1],&nindices);
//      free(biui0);
//      xi1indx+=ha_ndblocks[j1+begblock[rank]];
//    //}
//    //MPI_Barrier(PETSC_COMM_WORLD);
//  }
  free(xi1);
  free(insize);
  //free(insizeDD);
  ierr = PetscFree(submatC);
//  for(i=0; i<nmatint; i++)for(j=0;j<nreg;j++)MatDestroy(&submatBij[i*(nreg+1)+j][0]);
//  for(i=0; i<nmatin; i++)free(submatBij[i]);
//  free(submatBij);
  //for(i=0; i<nmatin; i++)MatDestroy(&submatCij[i]);
  //ierr = PetscFree(submatCij);
  //free(vnz);
  //free(begblock);
  printf("Solution calculation Rank %d time %f\n",rank,((double)(clock()-timestr))/CLOCKS_PER_SEC);
  //printf("Ok0 rank %d\n",rank);
  MPI_Barrier(PETSC_COMM_WORLD);
  //printf("Ok0 rank %d\n",rank);
  free(xd);
  timestr=clock();
  if(SORD==1)MPI_Allreduce(x0,x1,VecSize, MPI_DOUBLE, MPI_SUM,PETSC_COMM_WORLD);
  else MPI_Allreduce(x0,x1,VecSize, MPI_FLOAT, MPI_SUM,PETSC_COMM_WORLD);
  printf("Reduce solution Rank %d time %f\n",rank,((double)(clock()-timestr))/CLOCKS_PER_SEC);
  free(x0);
  //free(cirank);
  //free(birank);
  //free(aord);
  //free(ndblockinrank);
//  for (i=0; i<nmatin; i++) {
//    ierr = ISDestroy(&colindices[i]);
//    CHKERRQ(ierr);
//  }
//  PetscFree(colindices);
//  for (i=0; i<nmatin; i++) {
//    ierr = ISDestroy(&colindicesij[i]);
//    CHKERRQ(ierr);
//  }
//  PetscFree(colindicesij);
  printf("OK rui rank %d\n",rank);
  return 0;
}

int HaNDBBDParPreNos(Mat A, Vec b, ha_cgetype *x1, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, uvadd ndblock,uvadd nreg,uvadd ntime,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,uvdim laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol) {//,bool iter
  IS *rowindices,*colindices,*rowindicesbc,*colindicesbc1,*colindicesbc2,*Cindices,*Bindices,*BBindices,*BBiindices;//,*colindicesbcpm,*colindicesbcpm1
  const PetscInt *nindices;
  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nmatint,nmatinplust,nmatminust,nrowcolin,i,i1,j,j0,j1,j2,j3,j4,j5,j6,j7,l0,l1,l2,l3,l4,l5,rank1,proc1=0,nnzmax,j1nz,j1irnbs;//,sumrowcolin
  Mat *submatAij,*submatC,*submatD,submatBT,submatCT,*submatBB;//,*submatCij,*submatBij;,*submatB
  Vec vecxd;
  PetscInt *ai,*aj,*a1i,*a1j,*a2i,*a2j,*aic,*ajc;
  PetscReal cntl6in;
  uvadd lasize,ldsize;
  long int *biviindx1,*biviindx0,j7l;
  long int *obiviindx1= (long int*)calloc(1,sizeof(long int));
  long int *obiviindx0= (long int*)calloc(1,sizeof(long int));
  size_t freadresult,fwrt;
  obiviindx1[0]=-1;
  obiviindx0[0]=-1;
  PetscInt nrow,ncol,nz,nrowc,ncolc,nrowb,ncolb,nzc,nzv;
  MPI_Status   status;
  clock_t timestr,timeend;
  int j2int,la1;
  bool ifremove=false;
  PetscScalar *vals,*val1s,*val2s,*valsc,vecval;
  ha_cgetype vval;
  PetscErrorCode ierr;
  PetscViewer viewer;
  FILE* fp1,*fp2,*fp3,*presolfile;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  int *ndblockinrank= (int *) calloc (mpisize,sizeof(int));
  char filename[1024],rankname[1024],j1name[1024],tempchar[1024];//,fn1[1024],fn2[1024],fn3[1024]
  if(rank<10)strcpy(rankname,"000");
  if(rank<100&&rank>9)strcpy(rankname,"00");
  if(rank<1000&&rank>99)strcpy(rankname,"0");
  if(rank>=1000)rankname[0]='\0';
  sprintf(filename, "%d",rank);
  strcat(rankname,filename);
  int nmatinBB=1;
  nmatint=(uvadd)ntime/mpisize;
  nmatminust=nmatint;
  if(nmatint*mpisize<ntime)nmatinplust=nmatint+1;
  else nmatinplust=nmatint;
  for(i=0; i<mpisize; i++)if(rank+1<=ntime-mpisize*nmatint)nmatint++;
  nmatin=(nreg+1)*nmatint;
  nmatinplus=(nreg+1)*nmatinplust;
  begblock[rank]=nmatin;
  printf("rank %d nmatin %d nmatint %d nmplus %d\n",rank,nmatin,nmatint,nmatinplus);
  for(i=0; i<mpisize; i++) {
    j=rank;
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    j=j2;
    MPI_Bcast(&begblock[j2],1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  j2=begblock[0];
  begblock[0]=0;
  for(i=1; i<mpisize; i++) {
    j1=begblock[i];
    begblock[i]=begblock[i-1]+j2;
    j2=j1;
  }
  for(i=0; i<mpisize; i++)printf("rank %d beg block %d\n",rank,begblock[i]);
  j2=nmatin;
  MPI_Bcast(&j2,1, MPI_INT,0, PETSC_COMM_WORLD);
  ndblockinrank[0]=j2;
  for(i=1; i<mpisize; i++) {
    j2=nmatin;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    ndblockinrank[i]=ndblockinrank[i-1]+j2;
  }
  int *birank= (int *) calloc (nmatin,sizeof(int));
  int *cirank= (int *) calloc (nmatin,sizeof(int));
  //for(i=0; i<nmatin; i++)cirank[i]=-1;
  int *aord= (int *) calloc (nmatin,sizeof(int));
  for(i=0; i<ndblock; i++) {
    j2=i-ndblockinrank[rank];
    if(j2<0&&j2>=-nmatin) {
      j1=i/(nreg+1);
      j=i-(nreg+1)*j1;
      if(rank==0)j3=0;
      else j3=ndblockinrank[rank-1];
      if(j<nreg)cirank[j2+nmatin]=(nreg+1)*j1+nreg-j3;
      else cirank[j2+nmatin]=-1;
    }
  }
  j1=0;
  for(j=0; j<mpisize; j++)if(j<rank)j1+=ndblockinrank[j];
  for(i=0; i<nmatin; i++) {
    j2=(j1+i)/(nreg+1);
    if(i-j2*(nreg+1)!=0)aord[i]=1;
  }
  ierr = PetscMalloc(nmatin*sizeof(IS **),&rowindices);
  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*sizeof(IS **),&rowindicesbc);
//  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&colindices);
  CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbc1);
  //CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbc2);
  //CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*sizeof(IS **),&Cindices);
//  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatint*sizeof(IS **),&Bindices);
//  CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatinBB*sizeof(IS **),&BBindices);
//  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatAij);
  CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatint*sizeof(Mat *),&submatC);
  //CHKERRQ(ierr);
//  ierr = PetscMalloc(nmatinBB*sizeof(Mat *),&submatBB);
//  CHKERRQ(ierr);
  Mat** submatBij = (Mat**)calloc(nmatin,sizeof(Mat*));
  for (i=0; i<nmatin; i++) {
    submatBij[i] = (Mat*)calloc(1,sizeof(Mat));
  }
  Mat** submatCij = (Mat**)calloc(nmatin,sizeof(Mat*));
  for (i=0; i<nmatin; i++) {
    submatCij[i] = (Mat*)calloc(1,sizeof(Mat));
  }
//  Mat** submatB1 = (Mat**)calloc(nmatint,sizeof(Mat*));
//  for (i=0; i<nmatint; i++) {
//    submatB1[i] = (Mat*)calloc(1,sizeof(Mat));
//  }
//  Mat** submatB2 = (Mat**)calloc(nmatint,sizeof(Mat*));
//  for (i=0; i<nmatint; i++) {
//    submatB2[i] = (Mat*)calloc(1,sizeof(Mat));
//  }

//  ha_cgetype **yi2= (ha_cgetype**)calloc(nmatin,sizeof(ha_cgetype*));
//  for (i=0; i<nmatin; i++) {
//    yi2[i] = (ha_cgetype*)calloc(ha_ndblocks[i+begblock[rank]],sizeof(ha_cgetype));
//  }
//  sumrowcolin=0;
//  for(j=0; j<ndblock; j++)sumrowcolin+=ha_ndblocks[j];
//  PetscInt *offblock= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
//  PetscInt *offblockrow= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
//  offblock[0]=0;
//  for(j=0; j<ndblock; j++)offblock[j+1]=offblock[j]+countvarintra1[j+1]-countvarintra1[j]-ha_ndblocks[j];
//  offblockrow[0]=0;
//  for(j=0; j<ndblock; j++)offblockrow[j+1]=offblockrow[j]+counteqnoadd[j]-ha_ndblocks[j];

  char **fn01= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn01[i] = (char*)calloc(1024,sizeof(char));
  char **fn02= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn02[i] = (char*)calloc(1024,sizeof(char));
  char **fn03= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn03[i] = (char*)calloc(1024,sizeof(char));

  timestr=clock();
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      if(j1<10)strcpy(j1name,"000");
      if(j1<100&&j1>9)strcpy(j1name,"00");
      if(j1<1000&&j1>99)strcpy(j1name,"0");
      if(j1>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j1);
      strcat(j1name,filename);
      strcpy(filename,"_irnv");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn01[j1],filename);
      strcpy(filename,"_keep");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn02[j1],filename);
      strcpy(filename,"_vav");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn03[j1],filename);
    }
  }
  printf("OK012345 rank %d mpis %d sizeof int %ld size of long int %ld\n",rank,mpisize,sizeof(int),sizeof(long int));
//  for(j=0; j<mpisize; j++) {
//    //printf("j %d\n",j);
//    for(j3=0; j3<nmatinplus; j3++) {
//      if(j3>=nmatin)i=nmatin-1;
//      else i=j3;
//      //printf("rank %d OK1!!! i %d j %d nmatin %d\n",rank,i,j,nmatin);
//      bfirst=counteq[i+begblock[j]];
//      bend=ha_ndblocks[i+begblock[j]];
//      MPI_Bcast(&bfirst,1, MPI_INT,j, PETSC_COMM_WORLD);
//      MPI_Bcast(&bend,1, MPI_INT,j, PETSC_COMM_WORLD);
//      PetscInt *indices= (PetscInt *) calloc (bend,sizeof(PetscInt));
//      ha_cgetype *yi0 = (ha_cgetype*)calloc(bend,sizeof(ha_cgetype));
//      for(j1=0; j1<bend; j1++) { //ha_ndblocks[i+begblock[j]]
//        j2=ha_rows[bfirst+j1];//+bfirst;
//        if(j2>=Istart&&j2<Iend)indices[j1]=j2;
//        else indices[j1]=-1;
//      }
//      //printf("rank %d i %d j %d bf %d be %d\n",rank,i,j,bfirst,bend);
//      VecGetValues(b,bend,indices,yi0);
//      HaReduce(yi0,(forint)bend,mpisize,rank,j);
//      if(rank==j&&i<nmatin) {
//        memcpy (yi2[i],yi0,bend*sizeof(ha_cgetype));
//      }
//      free(indices);//1
//      free(yi0);//1
//    }
//  }
//  ha_cgetype **yi1= (ha_cgetype**)calloc(nmatint,sizeof(ha_cgetype*));
//  for (i=0; i<nmatint; i++) {
//    j1=0;
//    for (j=0; j<nreg+1; j++) {
//      j1+=ha_ndblocks[i*(nreg+1)+j+begblock[rank]];
//    }
//    yi1[i] = (ha_cgetype*)calloc(j1,sizeof(ha_cgetype));
//  }
//  for (i=0; i<nmatint; i++) {
//    j1=0;
//    for (j=0; j<nreg+1; j++) {
//      memcpy (&yi1[i][j1],yi2[i*(nreg+1)+j],ha_ndblocks[i*(nreg+1)+j+begblock[rank]]*sizeof(ha_cgetype));
//      free(yi2[i*(nreg+1)+j]);//1
//      j1+=ha_ndblocks[i*(nreg+1)+j+begblock[rank]];
//    }
//  }
//  free(yi2);//1
  printf("OK12345 rank %d\n",rank);
  j1=0;
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (ha_ndblocks[i+begblock[rank]],sizeof(PetscInt));
    //printf("ha_row 0 %d 588 %d\n",ha_rows[0],ha_rows[588]);
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_rows[bfirst+j];//+bfirst;
    }

    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
    bfirst=countvarintra1[i+begblock[rank]];
    //if(rank==0)printf("rank %d i %d bfirst %d block %d\n",rank,i,bfirst,ha_ndblocks[i+begblock[rank]]);
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_cols[bfirst+j];//+bfirst;
    }
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,colindices+i);
    free(indices);//1
  }
//  for(i=0; i<nmatint; i++) {
//    j2=0;
//    for(j1=0; j1<nreg+1; j1++)j2+=ha_ndblocks[i*(nreg+1)+j1+begblock[rank]];
//    PetscInt *indices= (PetscInt *) calloc (j2,sizeof(PetscInt));
//    j2=0;
//    for(j1=0; j1<nreg+1; j1++) {
//      bfirst=counteq[i*(nreg+1)+j1+begblock[rank]];
//      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
//        indices[j2]=ha_rows[bfirst+j];//+bfirst;
//        j2++;
//      }
//    }
//    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,rowindicesbc+i);
//    j2=0;
//    for(j1=0; j1<nreg; j1++) {
//      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
//      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
//        indices[j2]=ha_cols[bfirst+j];//+bfirst;
//        j2++;
//      }
//    }
//    //ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbc1+i);
//    j2=0;
//    for(j1=nreg; j1<nreg+1; j1++) {
//      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
//      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
//        indices[j2]=ha_cols[bfirst+j];//+bfirst;
//        j2++;
//      }
//    }
  //ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbc2+i);
//    free(indices);
//  }
//  printf("0okkkk\n");
//  PetscInt *indicesC= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
//  PetscInt *indicesB= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
//  for(i=0; i<ndblock; i++) {
//    bend=countvarintra1[i+1];
//    for(j=offblock[i]; j<offblock[i+1]; j++) {
//      j1=j-offblock[i];
//      indicesC[j]=ha_cols[bend-j1-1];//+bfirst;
//    }
//  }
//  j1=offblock[ndblock];
//  for(i=countvarintra1[ndblock]; i<VecSize; i++) {
//    indicesC[j1]=i;
//    j1++;
//  }
//  printf("1okkkk\n");
//  for(i=0; i<ndblock; i++) {
//    bend=counteq[i]+counteqnoadd[i];
//    for(j=offblockrow[i]; j<offblockrow[i+1]; j++) {
//      j1=j-offblockrow[i];
//      indicesB[j]=ha_rows[bend-j1-1];//+bfirst;
//    }
//  }
//  j1=offblockrow[ndblock];
//  for(i=counteq[ndblock-1]+counteqnoadd[ndblock-1]; i<VecSize; i++) {
//    indicesB[j1]=i;
//    j1++;
//  }
//  printf("2okkkk\n");
//  printf("rank1 %d\n",rank);
//  for(i=0; i<VecSize-sumrowcolin; i++)if(indicesB[i]>=Istart) {
//      j=i;
//      break;
//    }
//  for(j2=i; j2<VecSize-sumrowcolin; j2++)if(indicesB[j2]>=Iend) {
//      j1=j2;
//      break;
//    }
//  if(j2==VecSize-sumrowcolin)j1=VecSize-sumrowcolin;
//  printf("rank %d j %d j1 %d istart %d iend %d\n",rank,j,j1,Istart,Iend);
//  ha_cgetype *vecbiui= (ha_cgetype *) calloc (VecSize-sumrowcolin,sizeof(ha_cgetype));
//  VecGetValues(b,j1-j,&indicesB[j],&vecbiui[j]);//implicite yd save mem
  //ierr = VecDestroy(&b);
  //CHKERRQ(ierr);
//  printf("rank2 %d\n",rank);
//  for(i=0; i<nmatint; i++) {
//    ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesC,PETSC_COPY_VALUES,Cindices+i);
//    ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesB,PETSC_COPY_VALUES,Bindices+i);
//  }
//  free(indicesB);//1
//  free(offblock);//1
//  free(offblockrow);//1
//  PetscInt *indexBB= (PetscInt *) calloc (VecSize,sizeof(PetscInt));
//  for(i=0; i<VecSize; i++)indexBB[i]=i;
//  ISCreateGeneral(PETSC_COMM_SELF,VecSize,indexBB,PETSC_COPY_VALUES,BBindices);
  printf("Begin preparation rank %d\n",rank);
  ierr = MatCreateSubMatrices(A,nmatin,rowindices,colindices,MAT_INITIAL_MATRIX,&submatAij);
  //if(rank==0)printf("rank %d i %d bbrowij %d\n",rank,(nmatint-1)*(nreg+1)+nreg,submatAij[(nmatint-1)*(nreg+1)+nreg]->cmap->n);
  int insizes=17;
  int *insize=(int *) calloc (insizes*nmatin,sizeof(int));
  printf("rank %d A\n",rank);
  for(i=0; i<nmatminust; i++) {
    for(j=0; j<nreg; j++) {
      j1=i*(nreg+1)+j;
      j2=i*(nreg+1)+nreg;
      ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
      ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindices+j1,MAT_INITIAL_MATRIX,&submatBij[j1]);
      //printf("rank %d i %d j3 %d bbrowij %d\n",rank,i,j,submatCij[j1][0]->cmap->n);
    }
  }
  //printf("rank %d A1\n",rank);
  if(nmatint<nmatinplust) for(i=0; i<nreg; i++) {//to create again!
      MatDestroy(&submatCij[(nmatint-1)*(nreg+1)+i][0]);
      MatDestroy(&submatBij[(nmatint-1)*(nreg+1)+i][0]);
    }
  for(j=0; j<nreg; j++) {
    j1=(nmatint-1)*(nreg+1)+j;
    j2=(nmatint-1)*(nreg+1)+nreg;
    ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
    ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindices+j1,MAT_INITIAL_MATRIX,&submatBij[j1]);
    //printf("rank %d i %d j3 %d bbrowij %d\n",rank,i,j1,submatCij[j1][0]->cmap->n);
  }
  printf("rank %d A\n",rank);
  //ierr = MatCreateSubMatrices(A,nmatint,rowindicesbc,Cindices,MAT_INITIAL_MATRIX,&submatC);
  //printf("C\n");
//  ierr = MatCreateSubMatrices(A,nmatinBB,Bindices,BBindices,MAT_INITIAL_MATRIX,&submatBB);
  //ierr = MatDestroy(&A);
  //CHKERRQ(ierr);
//  PetscInt BBrow,BBcol;//,CCrow,CCcol;
//  MatGetSize(submatBB[0],&BBrow,&BBcol);
//  //MatGetSize(submatC[0],&CCrow,&CCcol);
//  printf("rank %d bbrow %d bbcol %d\n",rank,BBrow,BBcol);//,CCrow,CCcol);
//  PetscInt *indexBBi= (PetscInt *) calloc (BBrow,sizeof(PetscInt));
//  for(i=0; i<BBrow; i++)indexBBi[i]=i;
//  for(i=0; i<nmatint; i++) {
//    ierr = ISDestroy(&Bindices[i]);
//    ISCreateGeneral(PETSC_COMM_SELF,BBrow,indexBBi,PETSC_COPY_VALUES,Bindices+i);
//  }
  //ierr = MatCreateSubMatrices(submatBB[0],nmatint,Bindices,colindicesbc,MAT_INITIAL_MATRIX,&submatB);
  //char           ordering[256] = MATORDERINGRCM;
//  for(i=0; i<nmatint; i++) {
//    ierr = MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc1[i],MAT_INITIAL_MATRIX,&submatB1[i][0]);
//    ierr = ISDestroy(&colindicesbc1[i]);
//    ierr = MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc2[i],MAT_INITIAL_MATRIX,&submatB2[i][0]);
//    ierr = ISDestroy(&colindicesbc2[i]);
//  }
//  ierr = PetscFree(colindicesbc1);//1
//  ierr = PetscFree(colindicesbc2);//1
//  printf("B\n");
//  if(rank==mpisize-1)proc1=1;
//  ierr = PetscMalloc(proc1*sizeof(Mat *),&submatD);
//  CHKERRQ(ierr);
//  printf("D\n");
//  MatCreateSubMatrices(submatBB[0],proc1,Bindices,Cindices,MAT_INITIAL_MATRIX,&submatD);
  printf("End partitionning rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&rowindices[i]);//1
    CHKERRQ(ierr);
    ierr = ISDestroy(&colindices[i]);//1
    CHKERRQ(ierr);
  }
//  for (i=0; i<nmatint; i++) {
//    ierr = ISDestroy(&Cindices[i]);//1
//    CHKERRQ(ierr);
//    ierr = ISDestroy(&Bindices[i]);//1
//    CHKERRQ(ierr);
//  }
  PetscFree(colindices);//1
//  ierr = ISDestroy(&BBindices[0]);//1
  PetscFree(rowindices);//1

//  for (i=0; i<nmatint; i++) {
//    ierr = ISDestroy(&rowindicesbc[i]);//1
//    CHKERRQ(ierr);
//  }
//  PetscFree(rowindicesbc);//1

//  PetscFree(Cindices);//1
//  PetscFree(Bindices);//1
//  PetscFree(BBindices);//1
//  PetscFree(indexBB);//1
  //PetscFree(indexBBi);//1
//  ierr = MatDestroy(&submatBB[0]);//1
//  PetscFree(submatBB);//1
//  CHKERRQ(ierr);
//  printf("OK1 rank %d sumrowcolin %d\n",rank,sumrowcolin);
  //ha_cgetype *xi1 = (ha_cgetype*)calloc(sumrowcolin,sizeof(ha_cgetype));
  long int *bivinzrow=(long int *) calloc (1,sizeof(long int));
  PetscInt *bivinzcol=(PetscInt *) calloc (1,sizeof(PetscInt));
  PetscReal *ccolnorms=(PetscReal *) calloc (1,sizeof(PetscReal));
  //ccolnorms=realloc(ccolnorms,1*sizeof(PetscReal));
  int bivirowsize=1,bivicolsize=1,bbrowij,ddrowi;
  long int vecbivisize=0,halfvec;
  ha_cgetype *xi1point;
  long int *bivinzrow1=(long int *) calloc (1,sizeof(long int));
  int *bivinzcol1=(int *) calloc (1,sizeof(int));
  //ha_cgetype *vecbivi= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
  uvadd xi1indx=0;
  long int nz0,nz1,nz2,nz3=j,halfj2;
  long int lnz,lj2;

  for(j3=0; j3<nmatint; j3++) {
    bivirowsize=1;
    bivicolsize=1;
    bbrowij=submatBij[j3*(nreg+1)][0]->rmap->n;
    printf("j3 %d bbrowij %d\n",j3,bbrowij);
    ccolnorms=realloc(ccolnorms,bbrowij*sizeof(PetscReal));
    bivinzrow=realloc(bivinzrow,bbrowij*sizeof(long int));
    memset(bivinzrow,0,bbrowij*sizeof(long int));
    bivinzcol=realloc(bivinzcol,bbrowij*sizeof(PetscInt));
    memset(bivinzcol,0,bbrowij*sizeof(PetscInt));
    for(j1=0; j1<nreg; j1++) {
      MatGetColumnNorms(submatCij[j1+j3*(nreg+1)][0],NORM_1,ccolnorms);
      for(i=0; i<bbrowij; i++)if(ccolnorms[i]>0)bivinzcol[i]++;
      Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatBij[j1+j3*(nreg+1)][0]->data;//*aa=subA->data;
      ai= ab->i;
      nz=ab->nz;
      nrow=submatBij[j1+j3*(nreg+1)][0]->rmap->n;
      //printf("j1 %d nrow %d\n",j1,nrow);
      for(i=0; i<nrow-1; i++)if(ai[i]!=ai[i+1])bivinzrow[i]++;
      if(ai[i]<nz)bivinzrow[i]++;
    }
    ddrowi=submatBij[j3*(nreg+1)][0]->rmap->n;
    printf("ddrowi %d\n",ddrowi);
    for(i=0; i<ddrowi; i++) {
      if(bivinzrow[i]>0) {
        bivinzrow[i]=bivirowsize;
        bivirowsize++;
      }
      if(bivinzcol[i]>0) {
        bivinzcol[i]=bivicolsize;
        bivicolsize++;
      }
    }
    printf("j3 %d\n",j3);
    bivirowsize--;
    bivicolsize--;
    bivinzrow1=realloc (bivinzrow1,bivirowsize*sizeof(long int));
    bivinzcol1=realloc (bivinzcol1,bivicolsize*sizeof(int));
    j=0;
    j2=0;
    for(i=0; i<ddrowi; i++) {
      if(bivinzrow[i]>0) {
        bivinzrow1[j]=i*bbrowij;
        j++;
      }
      if(bivinzcol[i]>0) {
        bivinzcol1[j2]=i;
        j2++;
      }
      bivinzrow[i]=(bivinzrow[i]-1)*bivicolsize;//faster MSOL
    }
    vecbivisize=((long int)bivirowsize)*((long int)bivicolsize);
    printf("vecbivisize %ld bivirowsize %d bivicolsize %d\n",vecbivisize,bivirowsize,bivicolsize);
    if(vecbivisize<0)printf("Overflow detected, the interface problem is too big!!!!!!\n");
    //vecbivi=realloc (vecbivi,vecbivisize*sizeof(int));
    ha_cgetype *vecbivi= (ha_cgetype *) calloc (vecbivisize,sizeof(ha_cgetype));
    printf("vecbivisize %ld rank %d\n",vecbivisize,rank);
    xi1indx=0;
    for(j1=0; j1<nreg; j1++) {
      j4=j1+j3*(nreg+1);
      time(&timestr);//=clock();
      Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatAij[j4]->data;//*aa=subA->data;
      MatTranspose(submatCij[j4][0],MAT_INITIAL_MATRIX,&submatCT);
      MatDestroy(&submatCij[j4][0]);//1
      Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatCT->data;//*aa=subA->data;
      ai= aa->i;
      aj= aa->j;
      vals=aa->a;
      nz=aa->nz;
      nrow=submatAij[j4]->rmap->n;
      ncol=submatAij[j4]->cmap->n;
      lasize=ceil((laA/100.0)*nz);
      int *irn=(int *) calloc (lasize,sizeof(int));
      int *jcn=(int *) calloc (lasize,sizeof(int));
      ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
      insize[j4*insizes+13]=bivirowsize;
      insize[j4*insizes+14]=bivicolsize;
      insize[j4*insizes]=nrow;
      insize[j4*insizes+1]=ncol;
      insize[j4*insizes+2]=nz;
      aic= ac->i;
      ajc= ac->j;
      valsc=ac->a;
      nzc=ac->nz;
      nrowc=submatCT->rmap->n;
      ncolc=submatCT->cmap->n;
      insize[j4*insizes+3]=nrowc;
      insize[j4*insizes+4]=ncolc;
      memcpy (irn+nz,ai,(nrow+1)*sizeof(PetscInt));
      memcpy (jcn,aj,nz*sizeof(PetscInt));
      memcpy (values,vals,nz*sizeof(ha_cgetype));
      //printf("nz %d j %d\n",nz,j);
      insize[j4*insizes+2]=nz;
      insize[j4*insizes+5]=nzc;
      MatDestroy(&submatAij[j4]);//1
      insize[j4*insizes+9]=laA;
      insize[j4*insizes+10]=rank;
      insize[j4*insizes+11]=j4;
      MatGetSize(submatBij[j4][0],&nrowb,&ncolb);
      Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatBij[j4][0]->data;//*aa=subA->data;
      ai= ab->i;//irnbs+j1irnbs;
      aj= ab->j;//jcnbs+j1nz;
      vals=ab->a;//valbs+j1nz;
      nz=ab->nz;//nzbs[j1];
      insize[j4*insizes+6]=nrowb;//ncolb;
      insize[j4*insizes+7]=ncolb;//nrowb;
      insize[j4*insizes+8]=nz;
      prep48_msol_(insize+j4*insizes,irn,jcn,values,aic,ajc,valsc,ai,aj,vals,vecbivi,bivinzrow,bivinzcol);
      MatDestroy(&submatCT);
      MatDestroy(&submatBij[j4][0]);//1
      //free(yi1[j1]);
      free(irn);//1
      free(jcn);//1
      free(values);//1
      time(&timeend);
      //printf("t %d Submatrix %d rank %d calculation time %f\n",j3,j1,rank,difftime(timeend,timestr));
    }
    vecbivisize=((long int)bivirowsize)*((long int)bivicolsize);//(VecSize-sumrowcolin)*(VecSize-sumrowcolin);  //
    j=0;
    for(i=vecbivisize-1; i!=-1; --i)if(vecbivi[i]!=0)j++;
    //if(j==0)j=1;
    obiviindx1=realloc(obiviindx1,(j+1)*sizeof(long int));
    biviindx1=obiviindx1;
    biviindx1++;
    lj2=0;
    for(i=0; i<vecbivisize; ++i) {
      if(vecbivi[i]!=0) {
        vecbivi[lj2]=vecbivi[i];
        biviindx1[lj2]=bivinzcol1[i%bivicolsize]+bivinzrow1[i/bivicolsize];//[i-j1*bivicolsize]+bivinzrow1[j1];//*nrowb;
        ++lj2;
      }
    }
    bivinzrow1=realloc (bivinzrow1,1*sizeof(long int));
    bivinzcol1=realloc (bivinzcol1,1*sizeof(int));
    //if(j2==0){j2=1;biviindx1[0]=0;}
    vecbivi=realloc(vecbivi,lj2*sizeof(ha_cgetype));
    //printf("J2 %d\n",j2);
    //for(i=0;i<j2;i++)printf("aii %f indx %d\n",vecbivi[i],biviindx1[i]);
    nz0=lj2;
    j4=nreg+j3*(nreg+1);
    Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatAij[j4]->data;//*aa=subA->data;
    ai= aa->i;
    aj= aa->j;
    vals=aa->a;
    lnz=aa->nz;
    nrow=submatAij[j4]->rmap->n;
    ncol=submatAij[j4]->cmap->n;
    printf("D ncol %d nrow %d nz %ld\n",ncol,nrow,lnz);
    obiviindx0=realloc(obiviindx0,(lnz+1)*sizeof(long int));
    biviindx0=obiviindx0;
    biviindx0++;
    lj2=0;
    for(i=0; i<nrow-1; i++) {
      j7l=i*ncol;
      for(j=ai[i]; j<ai[i+1]; j++) {
        biviindx0[lj2]=aj[j]+j7l;
        lj2++;
      }
    }
    j7l=i*ncol;
    for(j=ai[i]; j<nz; j++) {
      biviindx0[lj2]=aj[j]+j7l;
      lj2++;
    }
    //printf("in C %ld\n",biviindx0[1]);
    my_spar_compl_(biviindx1,&nz0,biviindx0,&lnz,&nz1);
    printf("rank %d nz0 %ld nz %ld nz1 %ld\n",rank,nz0,lnz,nz1);
    vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
    int *irn1=(int *) calloc (nz1,sizeof(int));
    int *jcn1=(int *) calloc (nz1,sizeof(int));
    printf("nz1 %ld\n",nz1);
    //for(j=0;j<nz;j++)printf("rank %d j %d bivi %ld v %lf\n",rank,j,obiviindx0[j],vals[j]);
    my_spar_add4l_(vecbivi,obiviindx1,irn1,jcn1,&nz0,vals,obiviindx0,&lnz,&nz1,&ncol);
    printf("nz1 %ld\n",nz1);
    MatDestroy(&submatAij[j4]);//1
    printf("Arank %d nz0 %ld nz %ld nz1 %ld\n",rank,nz0,lnz,nz1);
    lnz=nz1;
    obiviindx1=realloc(obiviindx1,sizeof(long int));
    obiviindx0=realloc(obiviindx0,sizeof(long int));
    ldsize=ceil((laDi/100.0)*nz1);
    irn1=realloc(irn1,ldsize*sizeof(int));
    jcn1=realloc(jcn1,ldsize*sizeof(int));
    vecbivi=realloc(vecbivi,ldsize*sizeof(ha_cgetype));
    insize[j4*insizes]=nrow;
    insize[j4*insizes+1]=ncol;
    insize[j4*insizes+2]=lnz;
    insize[j4*insizes+10]=rank;
    insize[j4*insizes+11]=j4;
    insize[j4*insizes+9]=laDi;
    //for(i=0;i<nz;i++)printf("aii %f irn1 %d jcn %d\n",vecbivi[i],irn1[i],jcn[i]);
    //printf("ok!!!!!!!! vecbiuisize %d nrow %d\n",ncol,nrow);
    if(presol) {
      int *insized=(int *) calloc (5,sizeof(int));
      insized[0]=nrow;
      insized[1]=ncol;
      insized[2]=lnz;
      insized[4]=laDi;
      //if(rank==2)for(i=0;i<ncol;i++)printf("arank %d j1 %d col %d sct %d svb %d\n",rank,j3,jcn[i],(int)sizeof(PetscReal),(int)sizeof(ha_cgetype));
      if(cntl6==0&&SORD==0)cntl6in=0.3;
      else cntl6in=cntl6;
      spec51_rank_(insized,&cntl6in,irn1,jcn1,vecbivi);
      if(j3<10)strcpy(j1name,"000");
      if(j3<100&&j3>9)strcpy(j1name,"00");
      if(j3<1000&&j3>99)strcpy(j1name,"0");
      if(j3>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j3);
      strcat(j1name,filename);
      strcpy(filename,"_rank");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      printf("Rank %d RRRRR! ncol %d nrow %d rank %d nz %ld\n",rank,nrow,ncol,insized[3],lnz);
      fwrt=fwrite(insized, sizeof(int), 5, presolfile);
      if(fwrt!=5)printf("File write error.");
      fclose(presolfile);
      strcpy(filename,"_row");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fwrite(irn1, sizeof(int), nrow, presolfile);
      if( fwrt!= nrow) printf("File write error.");
      fclose(presolfile);
      strcpy(filename,"_col");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fwrite(jcn1, sizeof(int), ncol, presolfile);
      //if(rank==2)for(i=0;i<ncol;i++)printf("rank %d j1 %d col %d\n",rank,j3,jcn1[i]);
      if( fwrt!= ncol) printf("File write error.");
      fclose(presolfile);
      free(insized);//1
    }
    else {
      prep48_alu_(insize+j4*insizes,irn1,jcn1,vecbivi);
    }
    printf("ok!!!!!!!!\n");
    //vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
    free(vecbivi);//1
    free(irn1);//1
    free(jcn1);//1
  }
  //if(presol)return 0;
  free(ccolnorms);//1
  ierr = PetscFree(submatAij);//1
  printf("Rank %d End of preparation!\n",rank);
  MPI_Barrier(PETSC_COMM_WORLD);

  //bbrowij=submatB[0]->rmap->n;
  //printf("j3 %d bbrowij %d\n",j3,bbrowij);
  //ha_cgetype *vecbivi= (ha_cgetype *) calloc (vecbivisize,sizeof(ha_cgetype));
//  for(j1=0; j1<nmatint; j1++) {
//    MatDestroy(&submatCT);
//    free(yi1[j1]);//1
//    MatDestroy(&submatB1[j1][0]);//submatBT);//1
//    ierr = PetscFree(submatB1[j1]);//1
//    MatDestroy(&submatB2[j1][0]);//submatBT);//1
//    ierr = PetscFree(submatB2[j1]);//1
  //ierr = PetscFree(ai12);
  //ierr = PetscFree(aj12);
  //ierr = PetscFree(vals12);
//  }
  free(bivinzcol);//1
  free(bivinzrow);//1
//  ierr = PetscFree(submatB1);//1
//  ierr = PetscFree(submatB2);//1
//  CHKERRQ(ierr);
//  free(yi1);//1
  free(bivinzrow1);//1
  free(bivinzcol1);//1
  //vecbivi=realloc(vecbivi,j2*sizeof(ha_cgetype));
//  if(rank==mpisize-1) {
//    MatDestroy(&submatD[0]);//1
//    //free(irn1);
//    //free(jcn);
//    //free(insizeD);
//    //x0=(ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));//realloc (x0,*sizeof(ha_cgetype));
//  }
//  free(vecbiui);//1
  free(obiviindx1);//1
  free(obiviindx0);//1
  //free(vecbivi);//1
//  free(indicesC);//1
//  PetscFree(submatD);//1
//  for(j1=0; j1<nmatint; j1++) {
//    for(i=0; i<nreg; i++) {
//      MatDestroy(&submatCij[j1*(nreg+1)+i][0]);//1
//      MatDestroy(&submatBij[j1*(nreg+1)+i][0]);//1
//    }
//  }
  for(i=0; i<nmatin; i++) {
    PetscFree(submatCij[i]);//1
    PetscFree(submatBij[i]);//1
  }
//  ierr = PetscFree(colindicesbc1);//1
//  ierr = PetscFree(colindicesbc2);//1
  PetscFree(submatCij);//1
  PetscFree(submatBij);//1

//  free(xi1);//1
  free(insize);//1
  //ierr = PetscFree(submatC);//1
  free(begblock);//1
  timestr=clock();
  for(j1=0; j1<nmatin; j1++) {
    free(fn01[j1]);//1
    free(fn02[j1]);//1
    free(fn03[j1]);//1
  }
  free(fn01);//1
  free(fn02);//1
  free(fn03);//1
  free(cirank);//1
  free(birank);//1
  free(aord);//1
  free(ndblockinrank);//1
  return 0;
}


int HaNDBBDParSolNos(Mat A, Vec b, ha_cgetype *x1, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, uvadd ndblock,uvadd nreg,uvadd ntime,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,uvdim laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol) {//2 steps calculation, seperate presol. Need to run twice
  IS *rowindices,*colindices,*rowindicesbc,*colindicesbc1,*colindicesbc2,*Cindices,*Bindices,*BBindices,*BBiindices;//,*colindicesbcpm,*colindicesbcpm1
  const PetscInt *nindices;
  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nmatint,nmatinplust,nmatminust,nrowcolin,sumrowcolin,i,i1,j,j0,j1,j2,j3,j4,j5,j6,j7,l0,l1,l2,l3,l4,l5,rank1,proc1=0,nnzmax,j1nz,j1irnbs;
  Mat *submatAij,*submatC,*submatD,submatBT,submatCT,*submatBB;//,*submatCij,*submatBij;,*submatB
  Vec vecxd;
  PetscInt *ai,*aj,*a1i,*a1j,*a2i,*a2j,*aic,*ajc;
  PetscReal cntl6in;
  uvadd lasize,ldsize;
  long int *biviindx1,*biviindx0,j7l;
  long int *obiviindx1= (long int*)calloc(1,sizeof(long int));
  long int *obiviindx0= (long int*)calloc(1,sizeof(long int));
  size_t freadresult,fwrt;
  //printf("l l i %d l i %d i %d",sizeof(long long int),sizeof(long int),sizeof(int));
//  char tempchar[256];
  //obiviindx1=realloc(obiviindx1,sizeof(int));
  obiviindx1[0]=-1;
  //obiviindx0=realloc(obiviindx0,sizeof(int));
  obiviindx0[0]=-1;
  PetscInt nrow,ncol,nz,nrowc,ncolc,nrowb,ncolb,nzc,nzv;
  //forint iv,jv,la1;
  MPI_Status   status;
  clock_t timestr,timeend;
  //forint forintwrp;
  int j2int,la1;
  bool ifremove=false;
  //long long int nrow,ncol,nz,nrowc,ncolc,nzc,nzv;
  PetscScalar *vals,*val1s,*val2s,*valsc,vecval;
  ha_cgetype vval;
  PetscErrorCode ierr;
  PetscViewer viewer;
  FILE* fp1,*fp2,*fp3,*presolfile;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  int *ndblockinrank= (int *) calloc (mpisize,sizeof(int));
  char filename[1024],rankname[1024],j1name[1024],tempchar[1024];//,fn1[1024],fn2[1024],fn3[1024]
  if(rank<10)strcpy(rankname,"000");
  if(rank<100&&rank>9)strcpy(rankname,"00");
  if(rank<1000&&rank>99)strcpy(rankname,"0");
  if(rank>=1000)rankname[0]='\0';
  sprintf(filename, "%d",rank);
  strcat(rankname,filename);
  int nmatinBB=1;
  //reg matrix in one process only
  //nmatin=(uvadd)ndblock/mpisize;
  //if(nmatin*mpisize<ndblock)nmatinplus=nmatin+1;
  //else nmatinplus=nmatin;
  nmatint=(uvadd)ntime/mpisize;
  nmatminust=nmatint;
  if(nmatint*mpisize<ntime)nmatinplust=nmatint+1;
  else nmatinplust=nmatint;
  for(i=0; i<mpisize; i++)if(rank+1<=ntime-mpisize*nmatint)nmatint++;
  nmatin=(nreg+1)*nmatint;
  nmatinplus=(nreg+1)*nmatinplust;
  begblock[rank]=nmatin;
  printf("rank %d nmatin %d nmatint %d nmplus %d\n",rank,nmatin,nmatint,nmatinplus);
  for(i=0; i<mpisize; i++) {
    j=rank;
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    j=j2;
    MPI_Bcast(&begblock[j2],1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  j2=begblock[0];
  begblock[0]=0;
  for(i=1; i<mpisize; i++) {
    j1=begblock[i];
    begblock[i]=begblock[i-1]+j2;
    j2=j1;
  }
  for(i=0; i<mpisize; i++)printf("rank %d beg block %d\n",rank,begblock[i]);
  j2=nmatin;
  MPI_Bcast(&j2,1, MPI_INT,0, PETSC_COMM_WORLD);
  ndblockinrank[0]=j2;
  for(i=1; i<mpisize; i++) {
    j2=nmatin;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    ndblockinrank[i]=ndblockinrank[i-1]+j2;
  }
  int *birank= (int *) calloc (nmatin,sizeof(int));
  int *cirank= (int *) calloc (nmatin,sizeof(int));
  //for(i=0; i<nmatin; i++)cirank[i]=-1;
  int *aord= (int *) calloc (nmatin,sizeof(int));
  for(i=0; i<ndblock; i++) {
    j2=i-ndblockinrank[rank];
    if(j2<0&&j2>=-nmatin) {
      j1=i/(nreg+1);
      j=i-(nreg+1)*j1;
      if(rank==0)j3=0;
      else j3=ndblockinrank[rank-1];
      if(j<nreg)cirank[j2+nmatin]=(nreg+1)*j1+nreg-j3;
      else cirank[j2+nmatin]=-1;
    }
  }
  //printf("OK1!!!\n");
  j1=0;
  for(j=0; j<mpisize; j++)if(j<rank)j1+=ndblockinrank[j];
  for(i=0; i<nmatin; i++) {
    j2=(j1+i)/(nreg+1);
    if(i-j2*(nreg+1)!=0)aord[i]=1;
  }
  //printf("OK1a!!!\n");
  ierr = PetscMalloc(nmatin*sizeof(IS **),&rowindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&rowindicesbc);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&colindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbc1);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbc2);
  CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbcpm);
  //CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbcpm1);
  //CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&Cindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&Bindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatinBB*sizeof(IS **),&BBindices);
  CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatAij);
  //CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatint*sizeof(Mat *),&submatC);
  //CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatCij);
  //CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatint*sizeof(Mat *),&submatB);
  //CHKERRQ(ierr);
  //ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatBij);
  //CHKERRQ(ierr);
  //printf("OK2!!!\n");
  ierr = PetscMalloc(nmatinBB*sizeof(Mat *),&submatBB);
  CHKERRQ(ierr);
  Mat** submatBij = (Mat**)calloc(nmatin,sizeof(Mat*));
  for (i=0; i<nmatin; i++) {
    submatBij[i] = (Mat*)calloc(1,sizeof(Mat));
  }
  Mat** submatCij = (Mat**)calloc(nmatin,sizeof(Mat*));
  for (i=0; i<nmatin; i++) {
    submatCij[i] = (Mat*)calloc(1,sizeof(Mat));
  }
  Mat** submatB1 = (Mat**)calloc(nmatint,sizeof(Mat*));
  for (i=0; i<nmatint; i++) {
    submatB1[i] = (Mat*)calloc(1,sizeof(Mat));
  }
  Mat** submatB2 = (Mat**)calloc(nmatint,sizeof(Mat*));
  for (i=0; i<nmatint; i++) {
    submatB2[i] = (Mat*)calloc(1,sizeof(Mat));
  }
  //printf("rank %d OK3!!!\n",rank);

  ha_cgetype **yi2= (ha_cgetype**)calloc(nmatin,sizeof(ha_cgetype*));
  for (i=0; i<nmatin; i++) {
    //printf("i %d i+b %d s %d\n",i,i+begblock[rank],ha_ndblocks[i+begblock[rank]]);
    yi2[i] = (ha_cgetype*)calloc(ha_ndblocks[i+begblock[rank]],sizeof(ha_cgetype));
  }
  //printf("rank %d OK4!!!\n",rank);
  sumrowcolin=0;
  for(j=0; j<ndblock; j++)sumrowcolin+=ha_ndblocks[j];
  //printf("rank %d OK5!!!\n",rank);
  //PetscScalar *yd = (PetscScalar*)calloc(VecSize-sumrowcolin,sizeof(PetscScalar));
  PetscInt *offblock= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
  PetscInt *offblockrow= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
  //printf("rank %d OK6!!!\n",rank);
  offblock[0]=0;
  for(j=0; j<ndblock; j++)offblock[j+1]=offblock[j]+countvarintra1[j+1]-countvarintra1[j]-ha_ndblocks[j];
  //printf("rank %d OK7!!!\n",rank);
  offblockrow[0]=0;
  for(j=0; j<ndblock; j++)offblockrow[j+1]=offblockrow[j]+counteqnoadd[j]-ha_ndblocks[j];
  //printf("OK123\n");
  //printf("rank %d OK8!!!\n",rank);

  char **fn01= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn01[i] = (char*)calloc(1024,sizeof(char));
  char **fn02= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn02[i] = (char*)calloc(1024,sizeof(char));
  char **fn03= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn03[i] = (char*)calloc(1024,sizeof(char));
  //printf("rank %d OK9!!!\n",rank);

  timestr=clock();
  //printf("rank %d OK10!\n",rank);
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      if(j1<10)strcpy(j1name,"000");
      if(j1<100&&j1>9)strcpy(j1name,"00");
      if(j1<1000&&j1>99)strcpy(j1name,"0");
      if(j1>=1000)j1name[0]='\0';
      //printf("rank %d nmatin %d matinp %d j1 %d\n",rank,nmatin,nmatinplus,j1);
      sprintf(filename, "%d",j1);
      strcat(j1name,filename);
      strcpy(filename,"_irnv");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn01[j1],filename);
      strcpy(filename,"_keep");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn02[j1],filename);
      strcpy(filename,"_vav");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn03[j1],filename);
    }
  }
  //printf("rank %d OK11\n",rank);
  //int longi=3*pow(10,9);
  printf("OK012345 rank %d mpis %d sizeof int %ld size of long int %ld\n",rank,mpisize,sizeof(int),sizeof(long int));
  for(j=0; j<mpisize; j++) {
    //printf("j %d\n",j);
    for(j3=0; j3<nmatinplus; j3++) {
      if(j3>=nmatin)i=nmatin-1;
      else i=j3;
      //printf("rank %d OK1!!! i %d j %d nmatin %d\n",rank,i,j,nmatin);
      bfirst=counteq[i+begblock[j]];
      bend=ha_ndblocks[i+begblock[j]];
      MPI_Bcast(&bfirst,1, MPI_INT,j, PETSC_COMM_WORLD);
      MPI_Bcast(&bend,1, MPI_INT,j, PETSC_COMM_WORLD);
      PetscInt *indices= (PetscInt *) calloc (bend,sizeof(PetscInt));
      ha_cgetype *yi0 = (ha_cgetype*)calloc(bend,sizeof(ha_cgetype));
      for(j1=0; j1<bend; j1++) { //ha_ndblocks[i+begblock[j]]
        j2=ha_rows[bfirst+j1];//+bfirst;
        if(j2>=Istart&&j2<Iend)indices[j1]=j2;
        else indices[j1]=-1;
      }
      //printf("rank %d i %d j %d bf %d be %d\n",rank,i,j,bfirst,bend);
      VecGetValues(b,bend,indices,yi0);
      HaReduce(yi0,(forint)bend,mpisize,rank,j);
      if(rank==j&&i<nmatin) {
        memcpy (yi2[i],yi0,bend*sizeof(ha_cgetype));
      }
      free(indices);
      free(yi0);
    }
  }
  ha_cgetype **yi1= (ha_cgetype**)calloc(nmatint,sizeof(ha_cgetype*));
  for (i=0; i<nmatint; i++) {
    j1=0;
    for (j=0; j<nreg+1; j++) {
      j1+=ha_ndblocks[i*(nreg+1)+j+begblock[rank]];
    }
    yi1[i] = (ha_cgetype*)calloc(j1,sizeof(ha_cgetype));
  }
  for (i=0; i<nmatint; i++) {
    j1=0;
    for (j=0; j<nreg+1; j++) {
      memcpy (&yi1[i][j1],yi2[i*(nreg+1)+j],ha_ndblocks[i*(nreg+1)+j+begblock[rank]]*sizeof(ha_cgetype));
      free(yi2[i*(nreg+1)+j]);
      j1+=ha_ndblocks[i*(nreg+1)+j+begblock[rank]];
    }
  }
  free(yi2);
  printf("OK12345 rank %d\n",rank);
  j1=0;
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (ha_ndblocks[i+begblock[rank]],sizeof(PetscInt));
    //printf("ha_row 0 %d 588 %d\n",ha_rows[0],ha_rows[588]);
    //if(rank==0)printf("rank %d i %d bfirst %d block %d\n",rank,i,bfirst,ha_ndblocks[i+begblock[rank]]);
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_rows[bfirst+j];//+bfirst;
      //if(rank==0&&i==0)printf("indices %d rank %d\n",indices[j],ha_ndblocks[i+begblock[rank]]);
    }

//    bend=counteq[i+begblock[rank]]+counteqnoadd[i+begblock[rank]];
//    for(j=offblockrow[i+begblock[rank]]; j<offblockrow[i+1+begblock[rank]]; j++) {
//      j1=j-offblockrow[i+begblock[rank]];
//      indicesB[j]=ha_rows[bend-j1-1]+bfirst;
//    }
    //if(rank==0)printf("rank %d i %d ha_ndblocks[i+begblock[rank]] %d\n",rank,i,ha_ndblocks[i+begblock[rank]]);
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
    //if(cirank[i]>-1)ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
    bfirst=countvarintra1[i+begblock[rank]];
    //if(rank==0)printf("rank %d i %d bfirst %d block %d\n",rank,i,bfirst,ha_ndblocks[i+begblock[rank]]);
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_cols[bfirst+j];//+bfirst;
      //if(rank==0&&i==90)printf("indices %d rank %d\n",indices[j],ha_ndblocks[i+begblock[rank]]);
    }
//    bend=countvarintra1[i+1+begblock[rank]];
//    for(j=offblock[i+begblock[rank]]; j<offblock[i+1+begblock[rank]]; j++) {
//      j1=j-offblock[i+begblock[rank]];
//      indicesC[j]=ha_cols[bend-j1-1]+bfirst;
//      if(rank==1)printf("rank %d cindices %d bend %d j1 %d j %d\n",rank,indicesC[j],ha_cols[bend-j1-1],j1,j);
//    }
    //MPI_Bcast(indicesC[offblock[i+begblock[rank]]],(offblock[i+1+begblock[rank]]-offblock[i+begblock[rank]])*sizeof(PetscInt), MPI_BYTE,rank, PETSC_COMM_WORLD);
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,colindices+i);
    //if(cirank[i]>-1)ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,colindicesij+i);
    free(indices);
  }
  //MPI_Allreduce(MPI_IN_PLACE,indicesC,VecSize-sumrowcolin, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
  for(i=0; i<nmatint; i++) {
    j2=0;
    for(j1=0; j1<nreg+1; j1++)j2+=ha_ndblocks[i*(nreg+1)+j1+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (j2,sizeof(PetscInt));
    j2=0;
    //j3=0;
    for(j1=0; j1<nreg+1; j1++) {
      bfirst=counteq[i*(nreg+1)+j1+begblock[rank]];
      //printf("rank %d bfirst %d\n",rank,bfirst);
      //printf("ha_row 0 %d 588 %d\n",ha_rows[0],ha_rows[588]);
      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
        indices[j2]=ha_rows[bfirst+j];//+bfirst;
        //if(rank==0&&i==nmatint-1)printf("rank %d indx %d\n",rank,indices[j2]);
        j2++;
        //if(rank==0&&i==1)printf("indices %d rank %d\n",indices[j],ha_ndblocks[i+rank*nmatinplus]);
      }
    }
    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,rowindicesbc+i);
    //if(cirank[i]>-1)ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
//    j2=0;
//    for(j1=0; j1<nreg+1; j1++) {
//      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
//      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
//        indices[j2]=j2+bfirst;
//        j2++;
//      }
//    }
//    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbcpm+i);
//    j2=0;
//    for(j1=0; j1<nreg+1; j1++) {
//      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
//      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
//        indices[j2]=ha_cols[bfirst+j]-bfirst;//+bfirst;
//        //if(i==0&&j1==2)printf("j2 %d cindices %d cindx-1 %d\n",j2,indices[j2],indices[j2-1]);
//        //if(j2>0)if(indices[j2]<=indices[j2-1])printf("i %d j1 %d j2 %d indx %d\n",i,j1,j2,indices[j2]);
//        j2++;
//      }
//    }
//    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbcpm1+i);
    j2=0;
    for(j1=0; j1<nreg; j1++) {
      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
        indices[j2]=ha_cols[bfirst+j];//+bfirst;
        //if(i==0&&j1==2)printf("j2 %d cindices %d cindx-1 %d\n",j2,indices[j2],indices[j2-1]);
        //if(j2>0)if(indices[j2]<=indices[j2-1])printf("i %d j1 %d j2 %d indx %d\n",i,j1,j2,indices[j2]);
        j2++;
      }
    }
    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbc1+i);
    j2=0;
    for(j1=nreg; j1<nreg+1; j1++) {
      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
        indices[j2]=ha_cols[bfirst+j];//+bfirst;
        //if(i==0&&j1==2)printf("j2 %d cindices %d cindx-1 %d\n",j2,indices[j2],indices[j2-1]);
        //if(j2>0)if(indices[j2]<=indices[j2-1])printf("i %d j1 %d j2 %d indx %d\n",i,j1,j2,indices[j2]);
        j2++;
      }
    }
    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbc2+i);
    free(indices);
  }
  //printf("OK12345 rank %d\n",rank);
  //ISCopy(colindices,colindicesij);
  //ISCopy(rowindices,rowindicesij);
  printf("0okkkk\n");
  PetscInt *indicesC= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
  PetscInt *indicesB= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
  for(i=0; i<ndblock; i++) {
    bend=countvarintra1[i+1];
    //bfirst=countvarintra1[i];
    for(j=offblock[i]; j<offblock[i+1]; j++) {
      j1=j-offblock[i];
      indicesC[j]=ha_cols[bend-j1-1];//+bfirst;
      //if(i==98)printf("rank %d cindices %d bend %d j1 %d j %d\n",rank,indicesC[j],ha_cols[bend-j1-1],j1,j);
    }
  }
  j1=offblock[ndblock];
  //printf("rank %d j1 %d count %d vecs %d ndblock %d v-s %d\n",rank,j1,countvarintra1[ndblock],VecSize,ndblock,VecSize-sumrowcolin);
  for(i=countvarintra1[ndblock]; i<VecSize; i++) {
    indicesC[j1]=i;
    //printf("rank %d i %d indx %d\n",rank,i,indicesC[j1]);
    j1++;
  }
  printf("1okkkk\n");
  for(i=0; i<ndblock; i++) {
    //bfirst=counteq[i];
    bend=counteq[i]+counteqnoadd[i];
    for(j=offblockrow[i]; j<offblockrow[i+1]; j++) {
      j1=j-offblockrow[i];
      indicesB[j]=ha_rows[bend-j1-1];//+bfirst;
    }
  }
  j1=offblockrow[ndblock];
  for(i=counteq[ndblock-1]+counteqnoadd[ndblock-1]; i<VecSize; i++) {
    indicesB[j1]=i;
    j1++;
  }
  printf("2okkkk\n");
//  for(i=0; i<mpisize; i++) {
//    j=offblockrow[begblock[rank]];
//    j1=offblockrow[nmatin+begblock[rank]]-offblockrow[begblock[rank]];
//    MPI_Bcast(&j,1, MPI_INT,i, PETSC_COMM_WORLD);
//    MPI_Bcast(&j1,1, MPI_INT,i, PETSC_COMM_WORLD);
//    MPI_Bcast(&indicesB[j],j1, MPI_INT,i, PETSC_COMM_WORLD);
//  }
//  for(i=0; i<VecSize-sumrowcolin; i++)printf("rank %d indxb %d\n",rank,indicesB[i]);
  printf("rank1 %d\n",rank);
  for(i=0; i<VecSize-sumrowcolin; i++)if(indicesB[i]>=Istart) {
      j=i;
      break;
    }
  for(j2=i; j2<VecSize-sumrowcolin; j2++)if(indicesB[j2]>=Iend) {
      j1=j2;
      break;
    }
  if(j2==VecSize-sumrowcolin)j1=VecSize-sumrowcolin;
  printf("rank %d j %d j1 %d istart %d iend %d\n",rank,j,j1,Istart,Iend);
  ha_cgetype *vecbiui= (ha_cgetype *) calloc (VecSize-sumrowcolin,sizeof(ha_cgetype));
  VecGetValues(b,j1-j,&indicesB[j],&vecbiui[j]);//implicite yd save mem
  ierr = VecDestroy(&b);
  CHKERRQ(ierr);
  printf("rank2 %d\n",rank);
//  for(i=0; i<mpisize; i++) {
//    j=offblock[begblock[rank]];
//    j1=offblock[nmatin+begblock[rank]]-offblock[begblock[rank]];
//    MPI_Bcast(&j,1, MPI_INT,i, PETSC_COMM_WORLD);
//    MPI_Bcast(&j1,1, MPI_INT,i, PETSC_COMM_WORLD);
//    MPI_Bcast(&indicesC[j],j1, MPI_INT,i, PETSC_COMM_WORLD);
//  }
//  printf("OK132\n");
  for(i=0; i<nmatint; i++) {
    ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesC,PETSC_COPY_VALUES,Cindices+i);
    ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesB,PETSC_COPY_VALUES,Bindices+i);
  }
  free(indicesB);
  free(offblock);
  free(offblockrow);
  PetscInt *indexBB= (PetscInt *) calloc (VecSize,sizeof(PetscInt));
  for(i=0; i<VecSize; i++)indexBB[i]=i;
  ISCreateGeneral(PETSC_COMM_SELF,VecSize,indexBB,PETSC_COPY_VALUES,BBindices);
  printf("Begin preparation rank %d\n",rank);
  ierr = MatCreateSubMatrices(A,nmatin,rowindices,colindices,MAT_INITIAL_MATRIX,&submatAij);
  //if(rank==0)printf("rank %d i %d bbrowij %d\n",rank,(nmatint-1)*(nreg+1)+nreg,submatAij[(nmatint-1)*(nreg+1)+nreg]->cmap->n);
  int insizes=17;
  int *insize=(int *) calloc (insizes*nmatin,sizeof(int));
//  for(j1=0; j1<nmatin; j1++){
//    if(aord[j1]==1){
//      time(&timestr);//=clock();
//      Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatA[j1]->data;//*aa=subA->data;
//      ai= aa->i;
//      aj= aa->j;
//      vals=aa->a;
//      nz=aa->nz;
//      nrow=submatA[j1]->rmap->n;
//      ncol=submatA[j1]->cmap->n;
//      int *irn=(int *) calloc (laA*nz,sizeof(int));
//      int *jcn=(int *) calloc (laA*nz,sizeof(int));
//      ha_cgetype *values= (ha_cgetype *) calloc (laA*nz,sizeof(ha_cgetype));
//      insize[j1*insizes]=nrow;
//      insize[j1*insizes+1]=ncol;
//      insize[j1*insizes+2]=nz;
//      memcpy (irn+nz,ai,(nrow+1)*sizeof(PetscInt));
//      memcpy (jcn,aj,nz*sizeof(PetscInt));
//      memcpy (values,vals,nz*sizeof(ha_cgetype));
//      MatDestroy(&submatA[j1]);
//      insize[j1*insizes+9]=laA;
//      insize[j1*insizes+10]=rank;
//      insize[j1*insizes+11]=j1;
//      prep48_alu_(insize+j1*insizes,irn,jcn,values);
//      free(irn);
//      free(jcn);
//      free(values);
//      time(&timeend);
//      printf("Submatrix %d rank %d calculation time %f\n",j1,rank,difftime(timeend,timestr));
//    }
//  }

  printf("rank %d A\n",rank);
  for(i=0; i<nmatminust; i++) {
    for(j=0; j<nreg; j++) {
      j1=i*(nreg+1)+j;
      j2=i*(nreg+1)+nreg;
      ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
      ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindices+j1,MAT_INITIAL_MATRIX,&submatBij[j1]);
      //printf("rank %d i %d j3 %d bbrowij %d\n",rank,i,j,submatCij[j1][0]->cmap->n);
    }
  }
  //printf("rank %d A1\n",rank);
//  if(nmatin==nmatinplus)strcpy(tempchar,"MAT_INITIAL_MATRIX");
//  else strcpy(tempchar,"MAT_REUSE_MATRIX");
  if(nmatint<nmatinplust) for(i=0; i<nreg; i++) {//to create again
      MatDestroy(&submatCij[(nmatint-1)*(nreg+1)+i][0]);
      MatDestroy(&submatBij[(nmatint-1)*(nreg+1)+i][0]);
    }
  for(j=0; j<nreg; j++) {
    j1=(nmatint-1)*(nreg+1)+j;
    j2=(nmatint-1)*(nreg+1)+nreg;
    ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
    ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindices+j1,MAT_INITIAL_MATRIX,&submatBij[j1]);
    //printf("rank %d i %d j3 %d bbrowij %d\n",rank,i,j1,submatCij[j1][0]->cmap->n);
  }
  printf("rank %d A\n",rank);
  ierr = MatCreateSubMatrices(A,nmatint,rowindicesbc,Cindices,MAT_INITIAL_MATRIX,&submatC);
  printf("C\n");
  ierr = MatCreateSubMatrices(A,nmatinBB,Bindices,BBindices,MAT_INITIAL_MATRIX,&submatBB);
  ierr = MatDestroy(&A);
  CHKERRQ(ierr);
  PetscInt BBrow,BBcol;//,CCrow,CCcol;
  MatGetSize(submatBB[0],&BBrow,&BBcol);
  //MatGetSize(submatC[0],&CCrow,&CCcol);
  printf("rank %d bbrow %d bbcol %d\n",rank,BBrow,BBcol);//,CCrow,CCcol);
  PetscInt *indexBBi= (PetscInt *) calloc (BBrow,sizeof(PetscInt));
  for(i=0; i<BBrow; i++)indexBBi[i]=i;
  for(i=0; i<nmatint; i++) {
    ierr = ISDestroy(&Bindices[i]);
    ISCreateGeneral(PETSC_COMM_SELF,BBrow,indexBBi,PETSC_COPY_VALUES,Bindices+i);
  }
  //ierr = MatCreateSubMatrices(submatBB[0],nmatint,Bindices,colindicesbc,MAT_INITIAL_MATRIX,&submatB);
  //char           ordering[256] = MATORDERINGRCM;
  for(i=0; i<nmatint; i++) {
    ierr = MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc1[i],MAT_INITIAL_MATRIX,&submatB1[i][0]);
    ierr = MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc2[i],MAT_INITIAL_MATRIX,&submatB2[i][0]);
    //ierr = MatCreateSubMatrices(submatBB[0],1,Bindices+i,colindicesbc+i,MAT_INITIAL_MATRIX,&submatB[i]);
    //ierr = MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc[i],MAT_INITIAL_MATRIX,&submatB[i][0]);
    //MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbcpm[i],MAT_INITIAL_MATRIX,&submatBT);
    //MatGetOrdering(A,ordering,&Bindices[i],&colindicesbcpm1[i]);
    //ISSetPermutation(Bindices[i]);
    //ISSetPermutation(colindicesbcpm1[i]);
    //MatPermute(submatBT,Bindices[i],colindicesbcpm1[i],&submatB[i][0]);
    //MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc[i],&submatBT);
    //MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc[i],&submatB[i][0]);
    //MatGetLocalSubMatrix(submatBB[0],Bindices[i],colindicesbc[i],&submatBT);
    //MatGetLocalSubMatrix(submatBB[0],Bindices[i],colindicesbc[i],&submatB[i][0]);
    //MatCopy(submatBT,submatB[i][0],DIFFERENT_NONZERO_PATTERN);
    //MatConvert(submatBT,MATSEQAIJ,MAT_INITIAL_MATRIX,&submatB[i][0]);
    //MatDuplicate(submatBT,MAT_COPY_VALUES,&submatB[i][0]);
    //MatRestoreLocalSubMatrix(submatBB[0],Bindices[i],colindicesbc[i],&submatBT);
    //printf("iiiii %d\n",i);
    //ierr = ISDestroy(&colindicesbcpm[i]);
    //ierr = ISDestroy(&colindicesbcpm1[i]);
    //MatDestroy(&submatBT);
  }
  //PetscFree(colindicesbcpm);
  //PetscFree(colindicesbcpm1);
  printf("B\n");
  if(rank==mpisize-1)proc1=1;
  ierr = PetscMalloc(proc1*sizeof(Mat *),&submatD);
  CHKERRQ(ierr);
  printf("D\n");
  MatCreateSubMatrices(submatBB[0],proc1,Bindices,Cindices,MAT_INITIAL_MATRIX,&submatD);
  printf("End partitionning rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&rowindices[i]);
    CHKERRQ(ierr);
    ierr = ISDestroy(&colindices[i]);
    CHKERRQ(ierr);
    //ierr = ISDestroy(&rowindicesij[i]);
    //CHKERRQ(ierr);
    //ierr = ISDestroy(&colindicesij[i]);
    //CHKERRQ(ierr);
  }
  for (i=0; i<nmatint; i++) {
    ierr = ISDestroy(&Cindices[i]);
    CHKERRQ(ierr);
    ierr = ISDestroy(&Bindices[i]);
    CHKERRQ(ierr);
  }
  PetscFree(colindices);
  ierr = ISDestroy(&BBindices[0]);
  PetscFree(rowindices);
  //PetscFree(rowindicesij);
  //PetscFree(colindicesij);

  for (i=0; i<nmatint; i++) {
    ierr = ISDestroy(&rowindicesbc[i]);
    CHKERRQ(ierr);
    //ierr = ISDestroy(&colindicesbc[i]);
    //CHKERRQ(ierr);
  }
  PetscFree(rowindicesbc);
  //PetscFree(colindicesbc);

  PetscFree(Cindices);
  PetscFree(Bindices);
  PetscFree(BBindices);
  PetscFree(indexBB);
  PetscFree(indexBBi);
  ierr = MatDestroy(&submatBB[0]);
  PetscFree(submatBB);
  CHKERRQ(ierr);
  printf("OK1 rank %d sumrowcolin %d\n",rank,sumrowcolin);
  ha_cgetype *xi1 = (ha_cgetype*)calloc(sumrowcolin,sizeof(ha_cgetype));
  long int *bivinzrow=(long int *) calloc (1,sizeof(long int));
  PetscInt *bivinzcol=(PetscInt *) calloc (1,sizeof(PetscInt));
  PetscReal *ccolnorms=(PetscReal *) calloc (1,sizeof(PetscReal));
  //ccolnorms=realloc(ccolnorms,1*sizeof(PetscReal));
  int bivirowsize=1,bivicolsize=1,bbrowij,ddrowi;
  long int vecbivisize=0,li,lj;
  ha_cgetype *xi1point;
  long int *bivinzrow1=(long int *) calloc (1,sizeof(long int));
  long int *bivinzcol1=(long int *) calloc (1,sizeof(long int));
  //ha_cgetype *vecbivi= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
  uvadd xi1indx=0;
  long int nz0,nz1,nz2,nz3=j,halfj2;
  long int lj2,lnz;

  for(j3=0; j3<nmatint; j3++) {
    bivirowsize=1;
    bivicolsize=1;
    bbrowij=submatBij[j3*(nreg+1)][0]->rmap->n;
    printf("j3 %d bbrowij %d\n",j3,bbrowij);
    ccolnorms=realloc(ccolnorms,bbrowij*sizeof(PetscReal));
    bivinzrow=realloc(bivinzrow,bbrowij*sizeof(long int));
    memset(bivinzrow,0,bbrowij*sizeof(long int));
    bivinzcol=realloc(bivinzcol,bbrowij*sizeof(PetscInt));
    memset(bivinzcol,0,bbrowij*sizeof(PetscInt));
    for(j1=0; j1<nreg; j1++) {
      MatGetColumnNorms(submatCij[j1+j3*(nreg+1)][0],NORM_1,ccolnorms);
      for(i=0; i<bbrowij; i++)if(ccolnorms[i]>0)bivinzcol[i]++;
      Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatBij[j1+j3*(nreg+1)][0]->data;//*aa=subA->data;
      ai= ab->i;
      nz=ab->nz;
      nrow=submatBij[j1+j3*(nreg+1)][0]->rmap->n;
      //printf("j1 %d nrow %d\n",j1,nrow);
      for(i=0; i<nrow-1; i++)if(ai[i]!=ai[i+1])bivinzrow[i]++;
      if(ai[i]<nz)bivinzrow[i]++;
    }
    ddrowi=submatBij[j3*(nreg+1)][0]->rmap->n;
    printf("ddrowi %d\n",ddrowi);
    for(i=0; i<ddrowi; i++) {
      if(bivinzrow[i]>0) {
        bivinzrow[i]=bivirowsize;
        bivirowsize++;
      }
      if(bivinzcol[i]>0) {
        bivinzcol[i]=bivicolsize;
        bivicolsize++;
      }
    }
    printf("j3 %d\n",j3);
    bivirowsize--;
    bivicolsize--;
    bivinzrow1=realloc (bivinzrow1,bivirowsize*sizeof(long int));
    bivinzcol1=realloc (bivinzcol1,bivicolsize*sizeof(long int));
    j=0;
    j2=0;
    for(li=0; li<ddrowi; li++) {
      if(bivinzrow[li]>0) {
        bivinzrow1[j]=li*bbrowij;
        j++;
      }
      if(bivinzcol[li]>0) {
        bivinzcol1[j2]=li;
        j2++;
      }
      bivinzrow[li]=(bivinzrow[li]-1)*bivicolsize;//faster MSOL
    }
    vecbivisize=((long int)bivirowsize)*((long int)bivicolsize);
    printf("vecbivisize %ld bivirowsize %d bivicolsize %d\n",vecbivisize,bivirowsize,bivicolsize);
    if(vecbivisize<0)printf("Overflow detected, the interface problem is too big!!!!!!\n");
    //vecbivi=realloc (vecbivi,vecbivisize*sizeof(int));
    ha_cgetype *vecbivi= (ha_cgetype *) calloc (vecbivisize,sizeof(ha_cgetype));
    printf("vecbivisize %ld rank %d\n",vecbivisize,rank);
    xi1indx=0;
    for(j1=0; j1<nreg; j1++) {
      j4=j1+j3*(nreg+1);
      time(&timestr);//=clock();
      Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatAij[j4]->data;//*aa=subA->data;
      MatTranspose(submatCij[j4][0],MAT_INITIAL_MATRIX,&submatCT);
      Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatCT->data;//*aa=subA->data;
      ai= aa->i;
      aj= aa->j;
      vals=aa->a;
      nz=aa->nz;
      nrow=submatAij[j4]->rmap->n;
      ncol=submatAij[j4]->cmap->n;
      lasize=ceil((laA/100.0)*nz);
      int *irn=(int *) calloc (lasize,sizeof(int));
      int *jcn=(int *) calloc (lasize,sizeof(int));
      ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
      insize[j4*insizes+13]=bivirowsize;
      insize[j4*insizes+14]=bivicolsize;
      insize[j4*insizes]=nrow;
      insize[j4*insizes+1]=ncol;
      insize[j4*insizes+2]=nz;
      aic= ac->i;
      ajc= ac->j;
      valsc=ac->a;
      nzc=ac->nz;
      nrowc=submatCT->rmap->n;
      ncolc=submatCT->cmap->n;
      insize[j4*insizes+3]=nrowc;
      insize[j4*insizes+4]=ncolc;
      memcpy (irn+nz,ai,(nrow+1)*sizeof(PetscInt));
      memcpy (jcn,aj,nz*sizeof(PetscInt));
      memcpy (values,vals,nz*sizeof(ha_cgetype));
      //printf("nz %d j %d\n",nz,j);
      insize[j4*insizes+2]=nz;
      insize[j4*insizes+5]=nzc;
      MatDestroy(&submatAij[j4]);
      insize[j4*insizes+9]=laA;
      insize[j4*insizes+10]=rank;
      insize[j4*insizes+11]=j4;
      MatGetSize(submatBij[j4][0],&nrowb,&ncolb);
      Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatBij[j4][0]->data;//*aa=subA->data;
      ai= ab->i;//irnbs+j1irnbs;
      aj= ab->j;//jcnbs+j1nz;
      vals=ab->a;//valbs+j1nz;
      nz=ab->nz;//nzbs[j1];
      insize[j4*insizes+6]=nrowb;//ncolb;
      insize[j4*insizes+7]=ncolb;//nrowb;
      insize[j4*insizes+8]=nz;
      //for(i=0;i<nz;i++)++aj[i];
      //xi1point=xi1+xi1indx;
      //printf("rank %d j3 %d j1 %d\n",rank,j3,j1);
      prep48_msol_(insize+j4*insizes,irn,jcn,values,aic,ajc,valsc,ai,aj,vals,vecbivi,bivinzrow,bivinzcol);
      MatDestroy(&submatCT);
      //free(yi1[j1]);
      free(irn);
      free(jcn);
      free(values);
      time(&timeend);
      //printf("t %d Submatrix %d rank %d calculation time %f\n",j3,j1,rank,difftime(timeend,timestr));
    }
    vecbivisize=((long int)bivirowsize)*((long int)bivicolsize);//(VecSize-sumrowcolin)*(VecSize-sumrowcolin);  //
    lj=0;
    for(li=vecbivisize-1; li!=-1; --li)if(vecbivi[li]!=0)lj++;
    //if(j==0)j=1;
    obiviindx1=realloc(obiviindx1,(lj+1)*sizeof(long int));
    biviindx1=obiviindx1;
    biviindx1++;
    lj2=0;
    for(li=0; li<vecbivisize; ++li) {
      if(vecbivi[li]!=0) {
        vecbivi[lj2]=vecbivi[li];
        biviindx1[lj2]=bivinzcol1[li%bivicolsize]+bivinzrow1[li/bivicolsize];//[i-j1*bivicolsize]+bivinzrow1[j1];//*nrowb;
        ++lj2;
      }
    }
    bivinzrow1=realloc (bivinzrow1,1*sizeof(long int));
    bivinzcol1=realloc (bivinzcol1,1*sizeof(long int));
    //if(j2==0){j2=1;biviindx1[0]=0;}
    vecbivi=realloc(vecbivi,lj2*sizeof(ha_cgetype));
    //printf("J2 %d\n",j2);
    //for(i=0;i<j2;i++)printf("aii %f indx %d\n",vecbivi[i],biviindx1[i]);
    nz0=lj2;
    j4=nreg+j3*(nreg+1);
    Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatAij[j4]->data;//*aa=subA->data;
    ai= aa->i;
    aj= aa->j;
    vals=aa->a;
    lnz=aa->nz;
    nrow=submatAij[j4]->rmap->n;
    ncol=submatAij[j4]->cmap->n;
    printf("D ncol %d nrow %d nz %ld\n",ncol,nrow,lnz);
    obiviindx0=realloc(obiviindx0,(lnz+1)*sizeof(long int));
    biviindx0=obiviindx0;
    biviindx0++;
    lj2=0;
    for(li=0; li<nrow-1; li++) {
      j7l=li*ncol;
      for(j=ai[li]; j<ai[li+1]; j++) {
        biviindx0[lj2]=aj[j]+j7l;
        lj2++;
      }
    }
    j7l=li*ncol;
    for(j=ai[i]; j<nz; j++) {
      biviindx0[lj2]=aj[j]+j7l;
      lj2++;
    }
    //printf("in C %ld\n",biviindx0[1]);
    my_spar_compl_(biviindx1,&nz0,biviindx0,&lnz,&nz1);
    printf("rank %d nz0 %ld nz %ld nz1 %ld\n",rank,nz0,lnz,nz1);
    vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
    int *irn1=(int *) calloc (nz1,sizeof(int));
    int *jcn1=(int *) calloc (nz1,sizeof(int));
    printf("nz1 %ld\n",nz1);
    //for(j=0;j<nz;j++)printf("rank %d j %d bivi %ld v %lf\n",rank,j,obiviindx0[j],vals[j]);
    my_spar_add4l_(vecbivi,obiviindx1,irn1,jcn1,&nz0,vals,obiviindx0,&lnz,&nz1,&ncol);
    printf("nz1 %ld\n",nz1);
    MatDestroy(&submatAij[j4]);
    printf("Arank %d nz0 %ld nz %ld nz1 %ld\n",rank,nz0,lnz,nz1);
    lnz=nz1;
    obiviindx1=realloc(obiviindx1,sizeof(long int));
    obiviindx0=realloc(obiviindx0,sizeof(long int));
    ldsize=ceil((laDi/100.0)*nz1);
    irn1=realloc(irn1,ldsize*sizeof(int));
    jcn1=realloc(jcn1,ldsize*sizeof(int));
    vecbivi=realloc(vecbivi,ldsize*sizeof(ha_cgetype));
    memset(vecbivi+nz1,0,ceil((laDi-100)*nz1/100.0)*sizeof(ha_cgetype));
    memset(irn1+nz1,0,ceil((laDi-1)*nz1/100.0)*sizeof(int));
    insize[j4*insizes]=nrow;
    insize[j4*insizes+1]=ncol;
    insize[j4*insizes+2]=nz1;
    insize[j4*insizes+10]=rank;
    insize[j4*insizes+11]=j4;
    insize[j4*insizes+9]=laDi;
    //for(i=0;i<nz;i++)printf("aii %f irn1 %d jcn %d\n",vecbivi[i],irn1[i],jcn[i]);
    //printf("ok!!!!!!!! vecbiuisize %d nrow %d\n",ncol,nrow);
    printf("ok!!!!!!!! %d\n",presol);
    if(presol) {
      int *insized=(int *) calloc (5,sizeof(int));
      insized[0]=nrow;
      insized[1]=ncol;
      insized[2]=lnz;
      insized[4]=laDi;
      //if(rank==2)for(i=0;i<ncol;i++)printf("arank %d j1 %d col %d sct %d svb %d\n",rank,j3,jcn[i],(int)sizeof(PetscReal),(int)sizeof(ha_cgetype));
      if(cntl6==0&&SORD==0)cntl6in=0.3;
      else cntl6in=cntl6;
      spec51_rank_(insized,&cntl6in,irn1,jcn1,vecbivi);
      if(j3<10)strcpy(j1name,"000");
      if(j3<100&&j3>9)strcpy(j1name,"00");
      if(j3<1000&&j3>99)strcpy(j1name,"0");
      if(j3>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j3);
      strcat(j1name,filename);
      strcpy(filename,"_rank");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      printf("Rank %d RRRRR! ncol %d nrow %d rank %d nz %d\n",rank,nrow,ncol,insized[3],nz);
      fwrt=fwrite(insized, sizeof(int), 5, presolfile);
      if(fwrt!=5)printf("File write error.");
      fclose(presolfile);
      strcpy(filename,"_row");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fwrite(irn1, sizeof(int), nrow, presolfile);
      if( fwrt!= nrow) printf("File write error.");
      fclose(presolfile);
      strcpy(filename,"_col");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fwrite(jcn1, sizeof(int), ncol, presolfile);
      //if(rank==2)for(i=0;i<ncol;i++)printf("rank %d j1 %d col %d\n",rank,j3,jcn1[i]);
      if( fwrt!= ncol) printf("File write error.");
      fclose(presolfile);
      free(insized);
    }
    else {
      prep48_alu_(insize+j4*insizes,irn1,jcn1,vecbivi);
    }
    printf("ok!!!!!!!!\n");
    //vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
    free(vecbivi);
    free(irn1);
    free(jcn1);
  }
  if(presol)return 0;
  ierr = PetscFree(submatAij);
  printf("Rank %d End of preparation!\n",rank);
  MPI_Barrier(PETSC_COMM_WORLD);

  //bbrowij=submatB[0]->rmap->n;
  //printf("j3 %d bbrowij %d\n",j3,bbrowij);
  ccolnorms=realloc(ccolnorms,BBrow*sizeof(PetscReal));
  bivinzrow=realloc(bivinzrow,BBrow*sizeof(long int));
  memset(bivinzrow,0,BBrow*sizeof(long int));
  bivinzcol=realloc(bivinzcol,BBrow*sizeof(PetscInt));
  memset(bivinzcol,0,BBrow*sizeof(PetscInt));
  printf("rank %d OKKK BBrow %d\n",rank,BBrow);
  for(j1=0; j1<nmatint; j1++) {
    //if(j1<nmatint) {
    //printf("rank %d OKKK1 j1 %d\n",rank,j1);
    MatGetColumnNorms(submatC[j1],NORM_1,ccolnorms);
    //printf("OKKK1 j1 %d ccol %d\n",j1,submatC[j1]->cmap->n);
    //printf("j11111 %d\n",j1);
    for(i=0; i<BBrow; i++)if(ccolnorms[i]>0)bivinzcol[i]++;
    //printf("OKKK1 j1 %d\n",j1);
    Mat_SeqAIJ         *a1b=(Mat_SeqAIJ*)submatB1[j1][0]->data;//*aa=subA->data;
    ai= a1b->i;
    nz=a1b->nz;
    nrow=submatB1[j1][0]->rmap->n;
    for(i=0; i<nrow-1; i++)if(ai[i]!=ai[i+1])bivinzrow[i]++;
    if(ai[i]<nz)bivinzrow[i]++;

    Mat_SeqAIJ         *a2b=(Mat_SeqAIJ*)submatB2[j1][0]->data;//*aa=subA->data;
    ai= a2b->i;
    nz=a2b->nz;
    nrow=submatB2[j1][0]->rmap->n;
    for(j=0; j<nrow-1; j++)if(ai[j]!=ai[j+1])bivinzrow[j]++;
    if(ai[j]<nz)bivinzrow[j]++;
    nz=nz+a1b->nz;
    //nrow=nrow+submatB1[j1][0]->rmap->n;
    //printf("j11111 %d i %d\n",j,i);
    //}
  }
  printf("OKKK\n");
  free(ccolnorms);
  bivirowsize=0;
  bivicolsize=0;
  printf("BBrow %d vec-sum %ld\n",BBrow,VecSize-sumrowcolin);
  for(i=0; i<VecSize-sumrowcolin; i++) {
    if(bivinzrow[i]>0) {
      bivinzrow[i]=bivirowsize;
      bivirowsize++;
    }
    else bivinzrow[i]=-1;
    if(bivinzcol[i]>0) {
      bivinzcol[i]=bivicolsize;
      bivicolsize++;
    }
    else bivinzcol[i]=-1;
  }
  //bivirowsize--;
  //bivicolsize--;
  bivinzrow1=realloc (bivinzrow1,bivirowsize*sizeof(long int));//(int *) calloc (bivirowsize,sizeof(int));
  memset(bivinzrow1,0,bivirowsize*sizeof(long int));
  bivinzcol1=realloc (bivinzcol1,bivicolsize*sizeof(long int));//(int *) calloc (bivicolsize,sizeof(int));
  memset(bivinzcol1,0,bivicolsize*sizeof(long int));
  j=0;
  j2=0;
  long int lBBrow=BBrow;
  for(li=0; li<VecSize-sumrowcolin; li++) {
    if(bivinzrow[li]>=0) {
      bivinzrow1[j]=li*lBBrow;
      //printf("j %d bivinzrow1[j] %ld\n",j,bivinzrow1[j]);
      j++;
      bivinzrow[li]=bivinzrow[li]*bivicolsize;//faster MSOL
    }
    if(bivinzcol[li]>=0) {
      bivinzcol1[j2]=li;
      j2++;
    }
  }
  printf("j %d rows %d j2 %d cols %d\n",j,bivirowsize,j2,bivicolsize);
  vecbivisize=((long int)bivirowsize)*((long int)bivicolsize);
  if(vecbivisize<0)printf("Overflow detected, the interface problem is too big!!!!!!\n");
  ha_cgetype *vecbivi= (ha_cgetype *) calloc (vecbivisize,sizeof(ha_cgetype));
  //vecbivi=realloc(vecbivi,vecbivisize*sizeof(ha_cgetype));
  //memset(vecbivi,0,vecbivisize*sizeof(ha_cgetype));
  printf("vecbivisize %ld rank %d\n",vecbivisize,rank);
  //uvadd xi1indx=0;
  xi1indx=0;
  //j4=0;
  //int *insizeDD=(int *) calloc (insizes*nmatint,sizeof(int));
  int **irnereg= (int**)calloc(nreg+1,sizeof(int*));
  int **keepreg= (int**)calloc(nreg+1,sizeof(int*));
  ha_cgetype **valereg = (ha_cgetype**)calloc(nreg+1,sizeof(ha_cgetype*));
  //printf("OKt here11 nrowc %d ncolc %d!\n",nrowc,ncolc);
  for(i=0; i<nreg+1; i++) {
    irnereg[i]= (int*)calloc(1,sizeof(int));
    keepreg[i]= (int*)calloc(1,sizeof(int));
    valereg[i] = (ha_cgetype*)calloc(1,sizeof(ha_cgetype));
  }
  ha_cgetype cntl[10],rinfo[10],error1[3];
  int icntl[20],info[20],maxcolc1=0,maxcolc=0,maxrowc1=0;
  for(j1=0; j1<nmatint; j1++) {
    MatGetSize(submatC[j1],&maxrowc1,&maxcolc1);
    if(maxcolc<maxcolc1)maxcolc=maxcolc1;
    if(maxcolc<maxrowc1)maxcolc=maxrowc1;
  }
  ha_cgetype *w= (ha_cgetype*)calloc(4*maxcolc,sizeof(ha_cgetype));
  int *iw= (int*)calloc(maxcolc,sizeof(int));
  for(j1=0; j1<nmatint; j1++) {
    time(&timestr);//=clock();
    MatTranspose(submatC[j1],MAT_INITIAL_MATRIX,&submatCT);
    Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatCT->data;//*aa=subA->data;
    //insizeDD[j1*insizes+13]=bivirowsize;
    //insizeDD[j1*insizes+14]=bivicolsize;
    aic= ac->i;
    ajc= ac->j;
    valsc=ac->a;
    nzc=ac->nz;
    nrowc=submatCT->rmap->n;
    ncolc=submatCT->cmap->n;
    //insizeDD[j1*insizes+3]=nrowc;
    //insizeDD[j1*insizes+4]=ncolc;
    //insizeDD[j1*insizes+5]=nzc;
    //insizeDD[j1*insizes+10]=rank;
    //insizeDD[j1*insizes+11]=j1;
    //MatGetSize(submatB1[j1][0],&nrowb,&ncolb);
    Mat_SeqAIJ         *a1b=(Mat_SeqAIJ*)submatB1[j1][0]->data;//*aa=subA->data;
    a1i= a1b->i;//irnbs+j1irnbs;
    a1j= a1b->j;//jcnbs+j1nz;
    val1s=a1b->a;//valbs+j1nz;
    nz=a1b->nz;//nzbs[j1];
    //MatGetSize(submatB2[j1][0],&nrowb,&ncolb);
    Mat_SeqAIJ         *a2b=(Mat_SeqAIJ*)submatB2[j1][0]->data;//*aa=subA->data;
    a2i= a2b->i;//irnbs+j1irnbs;
    a2j= a2b->j;//jcnbs+j1nz;
    val2s=a2b->a;//valbs+j1nz;
    nz=nz+a2b->nz;//nzbs[j1];
    nrowb=submatB1[j1][0]->rmap->n;
    ncolb=submatB1[j1][0]->cmap->n;
    PetscInt *ai12= (PetscInt*)calloc(nrowb+1,sizeof(PetscInt));
    PetscInt *aj12= (PetscInt*)calloc(nz,sizeof(PetscInt));
    PetscScalar *vals12= (PetscScalar*)calloc(nz,sizeof(PetscScalar));
    for(i=0; i<nrowb-1; i++) {
      ai12[i]=a1i[i]+a2i[i];
      memcpy(&aj12[ai12[i]],&a1j[a1i[i]],(a1i[i+1]-a1i[i])*sizeof(PetscInt));
      memcpy(&vals12[ai12[i]],&val1s[a1i[i]],(a1i[i+1]-a1i[i])*sizeof(PetscScalar));
      j2=ai12[i]+a1i[i+1]-a1i[i];
      memcpy(&vals12[j2],&val2s[a2i[i]],(a2i[i+1]-a2i[i])*sizeof(PetscScalar));
      for(j=a2i[i]; j<a2i[i+1]; j++)aj12[j2+j-a2i[i]]=a2j[j]+ncolb;
    }
    ai12[i]=a1i[i]+a2i[i];
    ai12[i+1]=nz;
    nz=a1b->nz;
    memcpy(&aj12[ai12[i]],&a1j[a1i[i]],(nz-a1i[i])*sizeof(PetscInt));
    memcpy(&vals12[ai12[i]],&val1s[a1i[i]],(nz-a1i[i])*sizeof(PetscScalar));
    j2=ai12[i]+nz-a1i[i];
    nz=a2b->nz;
    memcpy(&vals12[j2],&val2s[a2i[i]],(nz-a2i[i])*sizeof(PetscScalar));
    for(j=a2i[i]; j<nz; j++)aj12[j2+j-a2i[i]]=a2j[j]+ncolb;
    nz=nz+a1b->nz;
    ncolb=ncolb+submatB2[j1][0]->cmap->n;
    //insizeDD[j1*insizes+6]=nrowb;//ncolb;
    //insizeDD[j1*insizes+7]=ncolb;//nrowb;
    //insizeDD[j1*insizes+8]=nz;
    //for(i=0;i<nz;i++)++aj[i];
    xi1point=xi1+xi1indx;
    printf("OK here j1 %d yi1 %ld xi1indx %ld a2j %d!\n",j1,sizeof(yi1[j1])/sizeof(PetscInt),xi1indx,ncolb);

    for(i=0; i<nreg+1; i++) {
      j2=j1*(nreg+1)+i;
      la1=ceil((insize[j2*insizes+9]/100.0)*insize[j2*insizes+2]);
      fp1 = fopen(fn01[j2], "rb");
      //printf("fn1 %s\n",fn01[j2]);
      if (fp1==NULL)printf("File opening error\n");
      irnereg[i] = realloc(irnereg[i],la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
      freadresult=fread(irnereg[i],sizeof(int),la1,fp1);
      fclose(fp1);
      fp2 = fopen(fn02[j2], "rb");
      if (fp2==NULL)printf("File opening error\n");
      keepreg[i] = realloc(keepreg[i],insize[j2*insizes+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
      freadresult=fread(keepreg[i],sizeof(int),insize[j2*insizes+12],fp2);
      fclose(fp2);
      fp3 = fopen(fn03[j2], "rb");
      if (fp3==NULL)printf("File opening error\n");
      valereg[i] = realloc(valereg[i],la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
      freadresult=fread(valereg[i],sizeof(ha_cgetype),la1,fp3);
      fclose(fp3);
    }
    printf("OK here11a nrowc %d ncolc %d rank %d!\n",nrowc,ncolc,rank);

    NDBBD_sol_nread(rank,j1*(nreg+1),nreg,insize,insizes,submatCij,submatBij,yi1[j1],xi1point,irnereg,keepreg,valereg,cntl,rinfo,error1,icntl,info,w,iw);
    //printf("Submatrix %d rank %d one mat calculation time %f\n",j1,rank,difftime(timeend,timestr));
    //if(j1==0)for(i=1;i<1000;i++)printf("j1 %d i %d xi %lf\n",j1,i,xi1point[i]);
    printf("OK here11b nrowc %d ncolc %d rank %d!\n",nrowc,ncolc,rank);
    ha_cgetype *bccol= (ha_cgetype*)calloc(ncolc,sizeof(ha_cgetype));
    for(i=0; i<nrowc-1; i++) {
      if(aic[i]<aic[i+1]) {
        memset(yi1[j1],0,ncolc*sizeof(ha_cgetype));
        for(j=aic[i]; j<aic[i+1]; j++) {
          yi1[j1][ajc[j]]=valsc[j];
        }
        //printf("OK here11 i %d nrowc %d!\n",i,nrowc);
        NDBBD_sol_nread(rank,j1*(nreg+1),nreg,insize,insizes,submatCij,submatBij,yi1[j1],bccol,irnereg,keepreg,valereg,cntl,rinfo,error1,icntl,info,w,iw);
        //printf("OK here112 i %d!\n",i);
        spar_vbiviadd_(bccol,bivinzcol+i,bivinzrow,&vecbivisize,&nrowb,&ncolc,&nz,ai12,aj12,vals12,vecbivi);
//        for(j2=0; j2<nrowb-1; j2++) {
//          for(j=ai[j2]; j<ai[j2+1]; j++) {
//            vecbivi[bivinzcol[i]+bivinzrow[j2]]-=vals[j]*bccol[aj[j]];//xi[i];
//            //printf("rank %d j1 %d j %d bccol %f",rank,j1,j,bccol[aj[j]]);
//          }
//        }
//        //printf("OK here112 i %d!\n",i);
//        for(j=ai[j2]; j<nz; j++) {
//          vecbivi[bivinzcol[i]+bivinzrow[j2]]-=vals[j]*bccol[aj[j]];//xi[i];
//          //printf("rank %d j1 %d j %d bccol %f",rank,j1,j,bccol[aj[j]]);
//        }
      }
    }
    printf("OK here11c nrowc %d ncolc %d rank %d!\n",nrowc,ncolc,rank);
    if(aic[i]<nzc) {
      memset(yi1[j1],0,ncolc*sizeof(ha_cgetype));
      for(j=aic[i]; j<nzc; j++) {
        yi1[j1][ajc[j]]=valsc[j];
      }
      NDBBD_sol_nread(rank,j1*(nreg+1),nreg,insize,insizes,submatCij,submatBij,yi1[j1],bccol,irnereg,keepreg,valereg,cntl,rinfo,error1,icntl,info,w,iw);
      spar_vbiviadd_(bccol,bivinzcol+i,bivinzrow,&vecbivisize,&nrowb,&ncolc,&nz,ai12,aj12,vals12,vecbivi);
//      for(j2=0; j2<nrowb-1; j2++) {
//        for(j=ai[j2]; j<ai[j2+1]; j++) {
//          vecbivi[bivinzcol[i]+bivinzrow[j2]]-=vals[j]*bccol[aj[j]];//xi[i];
//        }
//      }
//      for(j=ai[j2]; j<nz; j++) {
//        vecbivi[bivinzcol[i]+bivinzrow[j2]]-=vals[j]*bccol[aj[j]];//xi[i];
//      }
    }
    printf("OK here11d nrowc %d ncolc %d rank %d!\n",nrowc,ncolc,rank);
    free(bccol);
    //spec48_msol_(insize+j1*insizes,irn,jcn,values,yi1[j1],xi1point,aic,ajc,valsc,ai,aj,vals,vecbivi,bivinzrow,bivinzcol);
    MatDestroy(&submatCT);
    free(yi1[j1]);
    //Multiply Bi by ui:
    spar_mulmin_(xi1point,&nrowb,&nz,ai12,aj12,vals12,vecbiui);
//    for(i=0; i<nrowb-1; i++) {
//      for(j=ai[i]; j<ai[i+1]; j++) {
//        vecbiui[i]-=vals[j]*xi1point[aj[j]];//xi[i];
//      }
//    }
//    for(j=ai[i]; j<nz; j++) {
//      vecbiui[i]-=vals[j]*xi1point[aj[j]];//xi[i];
//      //printf("i %d ui %lf\n",aj[j],xi1point[aj[j]]);
//    }
    //if(j1==10)for(j=0; j<ncolc; j++)printf("rank %d j1 %d x1 %f\n",rank,j1,xi1point[j]);
    MatDestroy(&submatB1[j1][0]);//submatBT);
    ierr = PetscFree(submatB1[j1]);
    MatDestroy(&submatB2[j1][0]);//submatBT);
    ierr = PetscFree(submatB2[j1]);
    ierr = PetscFree(ai12);
    ierr = PetscFree(aj12);
    ierr = PetscFree(vals12);
    xi1indx+=ncolc;
//    for(i=0; i<nreg; i++) {
//      //MatDestroy(&submatBij[j1*(nreg+1)+i][0]);
//      MatDestroy(&submatCij[j1*(nreg+1)+i][0]);
//    }
    time(&timeend);
    printf("Submatrix %d rank %d calculation time %f\n",j1,rank,difftime(timeend,timestr));
    //j4+=ncolc;
  }
  free(w);
  free(iw);
  printf("Completed calculation of partitioned matrices! Rank %d\n",rank);
  for(i=0; i<nreg+1; i++) {
    free(irnereg[i]);//= (int*)calloc(1,sizeof(int));
    free(keepreg[i]);//= (int*)calloc(1,sizeof(int));
    free(valereg[i]);// = (ha_cgetype*)calloc(1,sizeof(ha_cgetype));
  }
  free(irnereg);//= (int**)calloc(nreg+1,sizeof(int*));
  free(keepreg);//= (int**)calloc(nreg+1,sizeof(int*));
  free(valereg);// = (ha_cgetype**)calloc(nreg+1,sizeof(ha_cgetype*));
  free(bivinzcol);
  free(bivinzrow);
  ierr = PetscFree(submatB1);
  ierr = PetscFree(submatB2);
  //ierr = PetscFree(submatBij);
//  for(i=0; i<nmatin; i++)free(submatCij[i]);
//  free(submatCij);
  CHKERRQ(ierr);
  free(yi1);
  //MPI_Barrier(PETSC_COMM_WORLD);
  printf("Completed calculation of partitioned matrices! Rank %d\n",rank);
  timestr=clock();

  j=0;
  for(i=vecbivisize-1; i!=-1; --i)if(vecbivi[i]!=0)j++;
  obiviindx1=realloc(obiviindx1,(j+1)*sizeof(long int));
  biviindx1=obiviindx1;
  biviindx1++;
  //int nz0,nz1,nz2,nz3=j,halfj2;
  lj2=0;
  printf("vbivi %ld colsize %d\n",vecbivisize,bivicolsize);
  for(i=0; i<vecbivisize; ++i) {
    if(vecbivi[i]!=0) {
      vecbivi[lj2]=vecbivi[i];
      biviindx1[lj2]=bivinzcol1[i%bivicolsize]+bivinzrow1[i/bivicolsize];//[i-j1*bivicolsize]+bivinzrow1[j1];//*nrowb;
      //if(rank==1)printf("j2 %d bivi %ld c %d r %ld\n",j2,biviindx1[j2],bivinzcol1[i%bivicolsize],bivinzrow1[i/bivicolsize]);
      //if(biviindx1[j2]>=vecbivisize)printf("b %d v %d\n",biviindx1[j2],vecbivisize);
      ++lj2;
    }
  }
  //for(i=0; i<bivirowsize; ++i)printf("OK here11! i %d j %d j2 %d\n",i,bivinzrow1[i],bivinzcol1[i]);
  free(bivinzrow1);
  free(bivinzcol1);
  printf("OK here11! j %d j2 %ld bivi %ld\n",j,lj2,biviindx1[lj2-1]);
  vecbivi=realloc(vecbivi,lj2*sizeof(ha_cgetype));
  printf("OK here11!\n");
  nz0=lj2;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len,name_len_max,name_beg,class_size,color,group_size,ha_id,group_size1,group_size11,ha_id1;
  MPI_Comm HA_COMM,HA1_COMM;
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
      if(strcmp(vec_pr_sname+i*name_len_max,vec_pr_name+j1*name_len_max)==0) {
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
    if(strcmp(processor_name,vec_pr_name+i*name_len_max)==0) {
      color=i;
      break;
    }
  }
  free(vec_pr_name);
  free(vec_pr_sname);
  MPI_Comm_split(PETSC_COMM_WORLD,color,rank,&HA_COMM);
  MPI_Comm_rank( HA_COMM, &ha_id);
  MPI_Comm_size(HA_COMM,&group_size);
  //printf("rank %d color %d group %d\n",rank,color,group_size);
  if(ha_id==group_size-1)color=1;
  else color=0;
  MPI_Comm_split(PETSC_COMM_WORLD,color,rank,&HA1_COMM);
  MPI_Comm_rank( HA1_COMM, &ha_id1);
  MPI_Comm_size(HA1_COMM,&group_size1);
  group_size11=group_size1;
  MPI_Bcast(&group_size11,1, MPI_INT,mpisize-1, PETSC_COMM_WORLD);
  //printf("rank %d color %d group11 %d\n",rank,color,group_size11);
  //printf("rank %d j2 %d!!!!!!!!!!!!\n",rank,j2);

  if(ha_id!=group_size-1) {
    //j2=nz3;
    MPI_Send(&lj2,1, MPI_LONG,group_size-1, 10, HA_COMM);
    if(SORD==1)MPI_Send(vecbivi,lj2, MPI_DOUBLE,group_size-1, 11, HA_COMM);
    else MPI_Send(vecbivi,lj2, MPI_FLOAT,group_size-1, 11, HA_COMM);
    MPI_Send(biviindx1,lj2, MPI_LONG,group_size-1, 12, HA_COMM);
    obiviindx1=realloc(obiviindx1,sizeof(long int));
    vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
  }
  if(ha_id==group_size-1) {
    for(j6=0; j6<group_size-1; j6++) {
      MPI_Recv(&lj2,1, MPI_LONG,  MPI_ANY_SOURCE,10, HA_COMM,&status);
      i=status.MPI_SOURCE;
      printf("rank1 %d j2 %ld!!!!!!!!!!!!\n",rank,lj2);
      obiviindx0=realloc(obiviindx0,(lj2+1)*sizeof(long int));
      biviindx0=obiviindx0;
      biviindx0++;
      printf("rank2 %d j2 %ld!!!!!!!!!!!!\n",rank,lj2);
      ha_cgetype *vecbivi0= (ha_cgetype *) calloc (lj2,sizeof(ha_cgetype));
      if(SORD==1)MPI_Recv(vecbivi0,lj2, MPI_DOUBLE,i,11, HA_COMM,&status);
      else MPI_Recv(vecbivi0,lj2, MPI_FLOAT,i,11, HA_COMM,&status);
      MPI_Recv(biviindx0,lj2, MPI_LONG,i,12, HA_COMM,&status);
      printf("rank3 %d j2 %ld!!!!!!!!!!!!\n",rank,lj2);
      my_spar_compl_(biviindx1,&nz0,biviindx0,&lj2,&nz1);
      printf("rank4 %d nz0 %ld j2 %ld nz1 %ld!!!!!!!!!!!!\n",rank,nz0,lj2,nz1);
      vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
      obiviindx1=realloc(obiviindx1,(nz1+1)*sizeof(long int));
      biviindx1=obiviindx1;
      biviindx1++;
      printf("rank5 %d j2 %ld!!!!!!!!!!!!\n",rank,lj2);
      my_spar_add3l_(vecbivi,obiviindx1,&nz0,vecbivi0,obiviindx0,&lj2,&nz1);
      printf("rank6 %d j2 %ld!!!!!!!!!!!!\n",rank,lj2);
      obiviindx0=realloc(obiviindx0,sizeof(long int));
      //free(obiviindx0);
      free(vecbivi0);
      printf("!!!!!!!!!!!! rank11 %d nz0 %ld j2 %ld nz1 %ld\n",rank,nz0,lj2,nz1);
      nz0=nz1;
    }
  }

  if(group_size11>1&&color==1) {
    nz3=nz0;
    if(ha_id1!=group_size1-1) {
      lj2=nz3;
      MPI_Send(&lj2,1, MPI_LONG,group_size1-1, 10, HA1_COMM);
      if(SORD==1)MPI_Send(vecbivi,lj2, MPI_DOUBLE,group_size1-1, 11, HA1_COMM);
      else MPI_Send(vecbivi,lj2, MPI_FLOAT,group_size1-1, 11, HA1_COMM);
      MPI_Send(biviindx1,lj2, MPI_LONG,group_size1-1, 12, HA1_COMM);
      obiviindx1=realloc(obiviindx1,sizeof(long int));
      vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
    }
    if(ha_id1==group_size1-1) {
      for(j6=0; j6<group_size1-1; j6++) {
        MPI_Recv(&lj2,1, MPI_LONG,  MPI_ANY_SOURCE,10, HA1_COMM,&status);
        i=status.MPI_SOURCE;
        obiviindx0=realloc(obiviindx0,(lj2+1)*sizeof(long int));
        biviindx0=obiviindx0;
        biviindx0++;
        ha_cgetype *vecbivi0= (ha_cgetype *) calloc (lj2,sizeof(ha_cgetype));
        if(SORD==1)MPI_Recv(vecbivi0,lj2, MPI_DOUBLE,i,11, HA1_COMM,&status);
        else MPI_Recv(vecbivi0,lj2, MPI_FLOAT,i,11, HA1_COMM,&status);
        MPI_Recv(biviindx0,lj2, MPI_LONG,i,12, HA1_COMM,&status);
        my_spar_compl_(biviindx1,&nz0,biviindx0,&lj2,&nz1);
        vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
        obiviindx1=realloc(obiviindx1,(nz1+1)*sizeof(long int));
        biviindx1=obiviindx1;
        biviindx1++;
        my_spar_add3l_(vecbivi,obiviindx1,&nz0,vecbivi0,obiviindx0,&lj2,&nz1);
        obiviindx0=realloc(obiviindx0,sizeof(long int));
        free(vecbivi0);
        printf("!!!!!!!!!!!! rank22 %d nz0 %ld j2 %ld nz1 %ld\n",rank,nz0,lj2,nz1);
        nz0=nz1;
      }
    }
  }
  forint vecbiuisize=0;
  vecbiuisize=VecSize-sumrowcolin;
  printf("after bivi rank %d\n",rank);
  HaReduceNoComp(vecbiui,vecbiuisize,mpisize,rank,mpisize-1);
  printf("Completed MPI_Reduce Operation! Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  if(rank!=mpisize-1)vecbiui=realloc(vecbiui,sizeof(ha_cgetype));//free(vecbiui);
  timestr=clock();
  ha_cgetype *xd;//= (ha_cgetype*)calloc(vecbiuisize,sizeof(ha_cgetype));
  ha_cgetype *x0;//,*vecbivi0;//= (ha_cgetype*)calloc(VecSize,sizeof(ha_cgetype));
  //int j7;
  //long int lnz;
  if(rank==mpisize-1) {
    Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatD[0]->data;//*aa=subA->data;
    ai= aa->i;
    aj= aa->j;
    vals=aa->a;
    lnz=aa->nz;
    nrow=submatD[0]->rmap->n;
    ncol=submatD[0]->cmap->n;
    printf("D ncol %d nrow %d nz %ld\n",ncol,nrow,lnz);
    obiviindx0=realloc(obiviindx0,(nz+1)*sizeof(long int));
    biviindx0=obiviindx0;
    biviindx0++;
    lj2=0;
    lBBrow=ncol;
    for(i=0; i<nrow-1; i++) {
      j7l=i*lBBrow;
      for(j=ai[i]; j<ai[i+1]; j++) {
        biviindx0[lj2]=aj[j]+j7l;
        lj2++;
      }
    }
    j7l=i*lBBrow;
    for(j=ai[i]; j<lnz; j++) {
      biviindx0[lj2]=aj[j]+j7l;
      lj2++;
    }
    my_spar_compl_(biviindx1,&nz0,biviindx0,&lnz,&nz1);
    printf("nz0 %ld nz %ld nz1 %ld\n",nz0,lnz,nz1);
    ldsize=ceil((laD/100.0)*nz1);
    vecbivi=realloc(vecbivi,ldsize*sizeof(ha_cgetype));
    int *irn1=(int *) calloc (nz1,sizeof(int));
    int *jcn=(int *) calloc (ldsize,sizeof(int));
    printf("nz1 %ld obin %ld ond1 %ld\n",nz1,obiviindx1[0],obiviindx0[0]);
    my_spar_add4l_(vecbivi,obiviindx1,irn1,jcn,&nz0,vals,obiviindx0,&lnz,&nz1,&ncol);
    printf("nz1 %ld\n",nz1);
    MatDestroy(&submatD[0]);
    lnz=nz1;
    obiviindx1=realloc(obiviindx1,sizeof(long int));
    obiviindx0=realloc(obiviindx0,sizeof(long int));
    irn1=realloc(irn1,ldsize*sizeof(int));
    int *insizeD=(int *) calloc (4,sizeof(int));
    insizeD[0]=nrow;
    insizeD[1]=ncol;
    insizeD[2]=lnz;
    insizeD[3]=laD;
    printf("ok\n");
    xd=(ha_cgetype *) calloc (vecbiuisize,sizeof(ha_cgetype));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
    printf("ok!!!!!!!! vecbiuisize %ld nrow %d\n",vecbiuisize,nrow);
    spec48_ssol2la_(insizeD,irn1,jcn,vecbivi,vecbiui,xd);
    printf("ok!!!!!!!!\n");

    vecbivi=realloc(vecbivi,1*sizeof(ha_cgetype));
    vecbiui=realloc(vecbiui,1*sizeof(ha_cgetype));
    //free(yd);
    free(irn1);
    free(jcn);
    //free(values);
    free(insizeD);
    x0=(ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));//realloc (x0,*sizeof(ha_cgetype));
    for(i=0; i<vecbiuisize; i++) {
      //VecGetValues(xsol,1,&i,&xd[i]);
      x0[indicesC[i]]=xd[i];
      //printf("i %d x1d %lf\n",i,xd[i]);
    }
    //VecView(xsol,0);
    //ierr = VecDestroy(&xsol);
    //CHKERRQ(ierr);
  }
  free(vecbiui);
  free(obiviindx1);
  free(obiviindx0);
  free(vecbivi);
  free(indicesC);
  if(rank!=mpisize-1) {
    xd=(ha_cgetype *) calloc (vecbiuisize,sizeof(ha_cgetype));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
    x0=(ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));//realloc (x0,*sizeof(ha_cgetype));
  }
  PetscFree(submatD);
  //PetscFree(submatCT);
  MPI_Barrier(PETSC_COMM_WORLD);
  if(SORD==1)MPI_Bcast(xd, vecbiuisize, MPI_DOUBLE,mpisize-1, PETSC_COMM_WORLD);
  else MPI_Bcast(xd, vecbiuisize, MPI_FLOAT,mpisize-1, PETSC_COMM_WORLD);
  printf("Interface problem solution Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  timestr=clock();
  xi1indx=0;


  for(j1=0; j1<nmatint; j1++) {
    time(&timestr);//=clock();
    Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatC[j1]->data;//*aa=subA->data;
    //printf("0 rank %d j1 %d vec %lf\n",rank,j1,vecbivi[22800]);
    ai= ac->i;
    aj= ac->j;
    vals=ac->a;
    nz=ac->nz;
    nrow=submatC[j1]->rmap->n;
    ncol=submatC[j1]->cmap->n;
    //printf("nrow %d\n",nrow);
    ha_cgetype *biui0= (ha_cgetype *) calloc (nrow,sizeof(ha_cgetype));
    ha_cgetype *be0 = (ha_cgetype*)calloc(nrow,sizeof(ha_cgetype));
    for(i=0; i<nrow-1; i++) {
      for(j=ai[i]; j<ai[i+1]; j++) {
        be0[i]+=vals[j]*xd[aj[j]];
      }
    }
    for(j=ai[i]; j<nz; j++) {
      be0[i]+=vals[j]*xd[aj[j]];
    }
    MatDestroy(&submatC[j1]);
    printf("OK here j1 %d yi1 %ld xi1indx %ld!\n",j1,sizeof(yi1[j1])/sizeof(PetscInt),xi1indx);
    ifremove=true;
    NDBBD_sol(rank,j1*(nreg+1),nreg,insize,insizes,submatCij,submatBij,be0,biui0,ifremove,fn01,fn02,fn03);
    printf("OK here11 nrowc %d ncolc %d!\n",nrow,ncol);
    free(be0);
    //MatDestroy(&submatB[j1]);//submatBT);
    //xi1indx+=nrow;
    for(i=0; i<nreg; i++) {
      MatDestroy(&submatCij[j1*(nreg+1)+i][0]);
      MatDestroy(&submatBij[j1*(nreg+1)+i][0]);
    }
    xi1point=xi1+xi1indx;
    ISGetIndices(colindicesbc1[j1],&nindices);
    ISGetSize(colindicesbc1[j1],&j);
    for(i=0; i<j; i++) {
      x0[nindices[i]]+=xi1point[i]-biui0[i];//vals[i]-biui0[i];//vecval;
    }
    ISRestoreIndices(colindicesbc1[j1],&nindices);
    ierr = ISDestroy(&colindicesbc1[j1]);
    ISGetIndices(colindicesbc2[j1],&nindices);
    for(j=i; j<nrow; j++) {
      x0[nindices[j-i]]+=xi1point[j]-biui0[j];//vals[i]-biui0[i];//vecval;
    }
    ISRestoreIndices(colindicesbc2[j1],&nindices);
    ierr = ISDestroy(&colindicesbc2[j1]);

    CHKERRQ(ierr);
    free(biui0);
    xi1indx+=nrow;
    time(&timeend);
    printf("Submatrix %d rank %d calculation time %f\n",j1,rank,difftime(timeend,timestr));
    //j4+=ncolc;
  }
  for(i=0; i<nmatin; i++) {
    PetscFree(submatCij[i]);
    PetscFree(submatBij[i]);
  }
  ierr = PetscFree(colindicesbc1);
  ierr = PetscFree(colindicesbc2);
  PetscFree(submatCij);
  PetscFree(submatBij);

//  for(j1=0; j1<nmatint; j1++) {
//    //if(j1<nmatin) {
//      ha_cgetype *biui0= (ha_cgetype *) calloc (ha_ndblocks[j1+begblock[rank]],sizeof(ha_cgetype));
//      fp1 = fopen(fn01[j1], "rb");
//      //printf("fn1 %s\n",fn01[j1]);
//      if (fp1==NULL)printf("File opening error\n");
//      fp2 = fopen(fn02[j1], "rb");
//      if (fp2==NULL)printf("File opening error\n");
//      fp3 = fopen(fn03[j1], "rb");
//      if (fp3==NULL)printf("File opening error\n");
//      la1=insize[j1*insizes+9]*insize[j1*insizes+2];
//      int *irne = (int*)calloc(la1,sizeof(int));
//      int *keep = (int*)calloc(insize[j1*insizes+12],sizeof(int));
//      ha_cgetype *vale = (ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
//      fread(irne,sizeof(int),la1,fp1);
//      fread(keep,sizeof(int),insize[j1*insizes+12],fp2);
//      fread(vale,sizeof(ha_cgetype),la1,fp3);
//      fclose(fp1);
//      fclose(fp2);
//      fclose(fp3);
//      remove(fn01[j1]);
//      remove(fn02[j1]);
//      remove(fn03[j1]);
//      //for(i=0;i<la;i++)printf("vav %f\n",vale[i]);
//      Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatC[j1]->data;//*aa=subA->data;
//      //printf("0 rank %d j1 %d vec %lf\n",rank,j1,vecbivi[22800]);
//      ai= ac->i;
//      aj= ac->j;
//      vals=ac->a;
//      nz=ac->nz;
//      nrow=submatC[j1]->rmap->n;
//      ncol=submatC[j1]->cmap->n;
//      //printf("nrow %d\n",nrow);
//      ha_cgetype *be0 = (ha_cgetype*)calloc(nrow,sizeof(ha_cgetype));
//      for(i=0; i<nrow-1; i++) {
//        for(j=ai[i]; j<ai[i+1]; j++) {
//          be0[i]+=vals[j]*xd[aj[j]];
//        }
//      }
//      for(j=ai[i]; j<nz; j++) {
//        be0[i]+=vals[j]*xd[aj[j]];
//      }
//      MatDestroy(&submatC[j1]);
//      spec48_esol_(insize+j1*insizes,irne,vale,keep,be0,biui0);
//      //printf("here!\n");
//      free(irne);
//      free(keep);
//      free(vale);
//      free(be0);
//      ISGetIndices(colindices[j1],&nindices);
//      xi1point=xi1+xi1indx;
//      for(i=0; i<ha_ndblocks[j1+begblock[rank]]; i++) {
//        x0[nindices[i]]+=xi1point[i]-biui0[i];//vals[i]-biui0[i];//vecval;
//      }
//      ISRestoreIndices(colindices[j1],&nindices);
//      free(biui0);
//      xi1indx+=ha_ndblocks[j1+begblock[rank]];
//    //}
//    //MPI_Barrier(PETSC_COMM_WORLD);
//  }
  free(xi1);
  free(insize);
  //free(insizeDD);
  ierr = PetscFree(submatC);
//  for(i=0; i<nmatint; i++)for(j=0;j<nreg;j++)MatDestroy(&submatBij[i*(nreg+1)+j][0]);
//  for(i=0; i<nmatin; i++)free(submatBij[i]);
//  free(submatBij);
  //for(i=0; i<nmatin; i++)MatDestroy(&submatCij[i]);
  //ierr = PetscFree(submatCij);
  //free(vnz);
  free(begblock);
  printf("Solution calculation Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  //printf("Ok0 rank %d\n",rank);
  MPI_Barrier(PETSC_COMM_WORLD);
  //printf("Ok0 rank %d\n",rank);
  free(xd);
  timestr=clock();
  if(SORD==1)MPI_Allreduce(x0,x1,VecSize, MPI_DOUBLE, MPI_SUM,PETSC_COMM_WORLD);
  else MPI_Allreduce(x0,x1,VecSize, MPI_FLOAT, MPI_SUM,PETSC_COMM_WORLD);
  printf("Reduce solution Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  free(x0);
  for(j1=0; j1<nmatin; j1++) {
    free(fn01[j1]);
    free(fn02[j1]);
    free(fn03[j1]);
  }
  free(fn01);
  free(fn02);
  free(fn03);
  free(cirank);
  free(birank);
  free(aord);
  free(ndblockinrank);
//  for (i=0; i<nmatin; i++) {
//    ierr = ISDestroy(&colindices[i]);
//    CHKERRQ(ierr);
//  }
//  PetscFree(colindices);
//  for (i=0; i<nmatin; i++) {
//    ierr = ISDestroy(&colindicesij[i]);
//    CHKERRQ(ierr);
//  }
//  PetscFree(colindicesij);
  printf("OK rui rank %d\n",rank);
  return 0;
}

bool NDBBD_sol(PetscInt rank, int begmat,int nreg,int * insize,int insizes, Mat **submatCij,Mat **submatBij,ha_cgetype *b,ha_cgetype *sol,bool ifremove,char** fn01,char** fn02, char** fn03) {
  FILE* fp1,*fp2,*fp3;
  PetscScalar *vals,*valsc,vecval;
  PetscInt *ai,*aj,*aic,*ajc,nrow,nz,maxrowcij;//,ncol
  PetscInt i,j,j1,j2,indx01,la1;
  ha_cgetype *b01,*b03,*sol1,*sol2;
  size_t freadresult;
//  char filename[1024],rankname[1024],j1name[1024],tempchar[1024];
  int *irne= NULL;//(int*)calloc(1,sizeof(int));
  int *keep= NULL;//(int*)calloc(1,sizeof(int));
  ha_cgetype *vale = NULL;//(ha_cgetype*)calloc(1,sizeof(ha_cgetype));
//  if(rank<10)strcpy(rankname,"000");
//  if(rank<100&&rank>9)strcpy(rankname,"00");
//  if(rank<1000&&rank>99)strcpy(rankname,"0");
//  sprintf(filename, "%d",rank);
//  strcat(rankname,filename);
  //char **fn01= (char**)calloc(nreg+1,sizeof(char*));
  //for (i=0; i<nreg+1; i++) fn01[i] = (char*)calloc(1024,sizeof(char));
  //char **fn02= (char**)calloc(nreg+1,sizeof(char*));
  //for (i=0; i<nreg+1; i++) fn02[i] = (char*)calloc(1024,sizeof(char));
  //char **fn03= (char**)calloc(nreg+1,sizeof(char*));
  //for (i=0; i<nreg+1; i++) fn03[i] = (char*)calloc(1024,sizeof(char));
  //indx01=0;
  maxrowcij=0;
  b01=b;
  sol1=sol;
  int fd1,fd2,fd3,frrsl1,frrsl2,frrsl3;
  for (j1=begmat; j1<nreg+begmat; j1++) {
    j2=j1*insizes;//-begmat;
    if((submatCij[j1][0]->rmap->n)>maxrowcij)maxrowcij=submatCij[j1][0]->rmap->n;
//    if(j1<10)strcpy(j1name,"000");
//    if(j1<100&&j1>9)strcpy(j1name,"00");
//    if(j1<1000&&j1>99)strcpy(j1name,"0");
//    sprintf(filename, "%d",j1);
//    strcat(j1name,filename);
//    strcpy(filename,"_irnv");
//    strcat(filename,rankname);
//    strcat(filename,j1name);
//    strcat(filename,".bin");
//    strcpy(fn01[j2],filename);
//    strcpy(filename,"_keep");
//    strcat(filename,rankname);
//    strcat(filename,j1name);
//    strcat(filename,".bin");
//    strcpy(fn02[j2],filename);
//    strcpy(filename,"_vav");
//    strcat(filename,rankname);
//    strcat(filename,j1name);
//    strcat(filename,".bin");
//    strcpy(fn03[j2],filename);
    #pragma omp critical(nsol)
    {
    if(isLinux==1) {
      la1=ceil((insize[j2+9]/100.0)*insize[j2+2]);
      irne = realloc(irne,la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
      keep = realloc(keep,insize[j2+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
      vale = realloc(vale,la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
      fd1 = open(fn01[j1],O_RDONLY);
      frrsl1=read(fd1,irne,sizeof(int)*la1);
      close(fd1);
      fd2 = open(fn02[j1],O_RDONLY);
      frrsl2=read(fd2,keep,insize[j2+12]*sizeof(int));
      close(fd2);
      fd3 = open(fn03[j1],O_RDONLY);
      frrsl3=read(fd3,vale,la1*sizeof(ha_cgetype));
      close(fd3);
    }
    else {
      la1=ceil((insize[j2+9]/100.0)*insize[j2+2]);
      irne = realloc(irne,la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
      keep = realloc(keep,insize[j2+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
      vale = realloc(vale,la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
      fp1 = fopen(fn01[j1], "rb");
      //printf("fn1 %s\n",fn01[j1]);
      if (fp1==NULL)printf("File opening error\n");
      freadresult=fread(irne,sizeof(int),la1,fp1);
      fclose(fp1);
      fp2 = fopen(fn02[j1], "rb");
      if (fp2==NULL)printf("File opening error\n");
      freadresult=fread(keep,sizeof(int),insize[j2+12],fp2);
      fclose(fp2);
      fp3 = fopen(fn03[j1], "rb");
      if (fp3==NULL)printf("File opening error\n");
      freadresult=fread(vale,sizeof(ha_cgetype),la1,fp3);
      fclose(fp3);
    }
    }
    //free(fn01[j1]);
    //free(fn02[j1]);
    //free(fn03[j1]);
    if(insize[j2+16]!=la1)insize[j2+16]=la1;
    spec48m_esol_(insize+j2,irne,vale,keep,b01,sol1);
    //indx01+=insize[j1*insizes];
    b01+=insize[j2];
    sol1+=insize[j2];
  }
  sol2=sol;
  for (j1=begmat; j1<nreg+begmat; j1++) {
    Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatBij[j1][0]->data;
    ai= ab->i;
    aj= ab->j;
    vals=ab->a;
    nz=ab->nz;
    nrow=submatBij[j1][0]->rmap->n;
    spar_mulmin_(sol2,&nrow,&nz,ai,aj,vals,b01);
//    for(i=0; i<nrow-1; i++) {
//      for(j=ai[i]; j<ai[i+1]; j++) {
//        b01[i]-=vals[j]*sol2[aj[j]];
//      }
//    }
//    for(j=ai[i]; j<nz; j++) {
//      b01[i]-=vals[j]*sol2[aj[j]];
//    }
    sol2+=insize[j1*insizes];
  }
  j2=j1;//-begmat;
//  if(j1<10)strcpy(j1name,"000");
//  if(j1<100&&j1>9)strcpy(j1name,"00");
//  if(j1<1000&&j1>99)strcpy(j1name,"0");
//  sprintf(filename, "%d",j1);
//  strcat(j1name,filename);
//  strcpy(filename,"_irnv");
//  strcat(filename,rankname);
//  strcat(filename,j1name);
//  strcat(filename,".bin");
//  strcpy(fn01[j2],filename);
//  strcpy(filename,"_keep");
//  strcat(filename,rankname);
//  strcat(filename,j1name);
//  strcat(filename,".bin");
//  strcpy(fn02[j2],filename);
//  strcpy(filename,"_vav");
//  strcat(filename,rankname);
//  strcat(filename,j1name);
//  strcat(filename,".bin");
//  strcpy(fn03[j2],filename);
    #pragma omp critical(nsol)
    {
  if(isLinux==1) {
    la1=ceil((insize[j1*insizes+9]/100.0)*insize[j1*insizes+2]);
    irne = realloc(irne,la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
    keep = realloc(keep,insize[j1*insizes+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
    vale = realloc(vale,la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
    fd1 = open(fn01[j2],O_RDONLY);
    frrsl1=read(fd1,irne,sizeof(int)*la1);
    close(fd1);
    fd2 = open(fn02[j2],O_RDONLY);
    frrsl2=read(fd2,keep,insize[j1*insizes+12]*sizeof(int));
    close(fd2);
    fd3 = open(fn03[j2],O_RDONLY);
    frrsl3=read(fd3,vale,la1*sizeof(ha_cgetype));
    close(fd3);
  }
  else {
    la1=ceil((insize[j1*insizes+9]/100.0)*insize[j1*insizes+2]);
    irne = realloc(irne,la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
    keep = realloc(keep,insize[j1*insizes+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
    vale = realloc(vale,la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
    fp1 = fopen(fn01[j2], "rb");
    //printf("fn1 %s\n",fn01[j1]);
    if (fp1==NULL)printf("File opening error\n");
    freadresult=fread(irne,sizeof(int),la1,fp1);
    fclose(fp1);
    fp2 = fopen(fn02[j2], "rb");
    if (fp2==NULL)printf("File opening error\n");
    freadresult=fread(keep,sizeof(int),insize[j1*insizes+12],fp2);
    fclose(fp2);
    fp3 = fopen(fn03[j2], "rb");
    if (fp3==NULL)printf("File opening error\n");
    freadresult=fread(vale,sizeof(ha_cgetype),la1,fp3);
    fclose(fp3);
  }
    }
  if(ifremove)remove(fn01[j2]);
  if(ifremove)remove(fn02[j2]);
  if(ifremove)remove(fn03[j2]);
//  free(fn01[j2]);
//  free(fn02[j2]);
//  free(fn03[j2]);
  if(insize[j1*insizes+16]!=la1)insize[j1*insizes+16]=la1;
  spec48m_esol_(insize+j1*insizes,irne,vale,keep,b01,sol1);
  sol2=sol;
  b03=b;
  ha_cgetype *b02 = (ha_cgetype*)calloc(maxrowcij,sizeof(ha_cgetype));
  for (j1=begmat; j1<nreg+begmat; j1++) {
    j2=j1;//-begmat;
    #pragma omp critical(nsol)
    {
    if(isLinux==1) {
      la1=ceil((insize[j1*insizes+9]/100.0)*insize[j1*insizes+2]);
      irne = realloc(irne,la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
      keep = realloc(keep,insize[j1*insizes+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
      vale = realloc(vale,la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
      fd1 = open(fn01[j2],O_RDONLY);
      frrsl1=read(fd1,irne,sizeof(int)*la1);
      close(fd1);
      fd2 = open(fn02[j2],O_RDONLY);
      frrsl2=read(fd2,keep,insize[j1*insizes+12]*sizeof(int));
      close(fd2);
      fd3 = open(fn03[j2],O_RDONLY);
      frrsl3=read(fd3,vale,la1*sizeof(ha_cgetype));
      close(fd3);
    }
    else {
      la1=ceil((insize[j1*insizes+9]/100.0)*insize[j1*insizes+2]);
      irne = realloc(irne,la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
      keep = realloc(keep,insize[j1*insizes+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
      vale = realloc(vale,la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
      fp1 = fopen(fn01[j2], "rb");
      if (fp1==NULL)printf("File opening error\n");
      freadresult=fread(irne,sizeof(int),la1,fp1);
      fclose(fp1);
      fp2 = fopen(fn02[j2], "rb");
      if (fp2==NULL)printf("File opening error\n");
      freadresult=fread(keep,sizeof(int),insize[j1*insizes+12],fp2);
      fclose(fp2);
      fp3 = fopen(fn03[j2], "rb");
      if (fp3==NULL)printf("File opening error\n");
      freadresult=fread(vale,sizeof(ha_cgetype),la1,fp3);
      fclose(fp3);
    }
    }
    if(ifremove)remove(fn01[j2]);
    if(ifremove)remove(fn02[j2]);
    if(ifremove)remove(fn03[j2]);
//    free(fn01[j2]);
//    free(fn02[j2]);
//    free(fn03[j2]);
    Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatCij[j1][0]->data;//*aa=subA->data;
    //printf("0 rank %d j1 %d vec %lf\n",rank,j1,vecbivi[22800]);
    ai= ac->i;
    aj= ac->j;
    vals=ac->a;
    nz=ac->nz;
    nrow=submatCij[j1][0]->rmap->n;
    //ncol=submatCij[j1][0]->cmap->n;
    //memset(b01,0,insize[j1*insizes+3]*sizeof(ha_cgetype));
//    ha_cgetype *b02 = (ha_cgetype*)calloc(nrow,sizeof(ha_cgetype));
    //ha_cgetype *b03 = (ha_cgetype*)calloc(nrow,sizeof(ha_cgetype));
    spar_mulnoadd_(sol1,&nrow,&nz,ai,aj,vals,b02);
//    for(i=0; i<nrow-1; i++) {
//      //b01[i]=0;
//      for(j=ai[i]; j<ai[i+1]; j++) {
//        b02[i]+=vals[j]*sol1[aj[j]];
//      }
//    }
//    //b01[i]=0;
//    for(j=ai[i]; j<nz; j++) {
//      b02[i]+=vals[j]*sol1[aj[j]];
//    }
    //MatDestroy(&submatCij[j1]);
    if(insize[j1*insizes+16]!=la1)insize[j1*insizes+16]=la1;
    spec48m_esol_(insize+j1*insizes,irne,vale,keep,b02,b03);
//    free(b02);
    for(j=0; j<nrow; j++)sol2[j]-=b03[j];
    //free(b03);
    sol2+=nrow;
    //printf("nrow %d crow %d ccol %d arow %d\n",nrow,insize[j1*insizes+3],insize[j1*insizes+4],insize[j1*insizes]);
    //b01+=insize[j1*insizes+3];
  }
  free(b02);
  free(irne);
  free(keep);
  free(vale);
//  free(fn01);
//  free(fn02);
//  free(fn03);
  //printf("HHere\n");
  return true;
}

bool NDBBD_sol_nread(PetscInt rank, int begmat,int nreg,int * insize,int insizes, Mat **submatCij,Mat **submatBij,ha_cgetype *b,ha_cgetype *sol,int** irnereg,int** keepreg,ha_cgetype** valereg,ha_cgetype *cntl,ha_cgetype *rinfo,ha_cgetype *error1,int *icntl,int *info,ha_cgetype *w,int *iw) {
  //FILE* fp1,*fp2,*fp3;
  PetscScalar *vals,*valsc;//,vecval;
  PetscInt *ai,*aj,*aic,*ajc,nrow,nz;//,ncol
  PetscInt i,j,j1,j2;//,indx01,la1;
  ha_cgetype *b01,*b03,*sol1,*sol2;//,*b02
  b01=b;
  sol1=sol;
  //printf("OKKK111222\n");
  for (j1=begmat; j1<nreg+begmat; j1++) {
    j2=j1-begmat;//-begmat;
    j=j1*insizes;
    //printf("j2 %d\n",j2);
    spec48_rpesol_(insize+j,irnereg[j2],valereg[j2],keepreg[j2],b01,sol1,cntl,rinfo,error1,icntl,info,w,iw);
    b01+=insize[j];
    sol1+=insize[j];
  }
  //printf("OKKK111222\n");
  sol2=sol;
  for (j1=begmat; j1<nreg+begmat; j1++) {
    Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatBij[j1][0]->data;
    ai= ab->i;
    aj= ab->j;
    vals=ab->a;
    nz=ab->nz;
    nrow=submatBij[j1][0]->rmap->n;
    spar_mulmin_(sol2,&nrow,&nz,ai,aj,vals,b01);
//    for(i=0; i<nrow-1; i++) {
//      for(j=ai[i]; j<ai[i+1]; j++) {
//        b01[i]-=vals[j]*sol2[aj[j]];
//      }
//    }
//    for(j=ai[i]; j<nz; j++) {
//      b01[i]-=vals[j]*sol2[aj[j]];
//    }
    sol2+=insize[j1*insizes];
  }
  j2=j1-begmat;//-begmat;
  spec48_rpesol_(insize+j1*insizes,irnereg[j2],valereg[j2],keepreg[j2],b01,sol1,cntl,rinfo,error1,icntl,info,w,iw);
  sol2=sol;
  //printf("OKKK111222\n");
  b03=b;
  for (j1=begmat; j1<nreg+begmat; j1++) {
    j2=j1-begmat;//-begmat;
    Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatCij[j1][0]->data;//*aa=subA->data;
    ai= ac->i;
    aj= ac->j;
    vals=ac->a;
    nz=ac->nz;
    nrow=submatCij[j1][0]->rmap->n;
    ha_cgetype *b02 = (ha_cgetype*)calloc(nrow,sizeof(ha_cgetype));
    //ha_cgetype *b03 = (ha_cgetype*)calloc(nrow,sizeof(ha_cgetype));
    spar_muladd_(sol1,&nrow,&nz,ai,aj,vals,b02);
//    for(i=0; i<nrow-1; i++) {
//      for(j=ai[i]; j<ai[i+1]; j++) {
//        b02[i]+=vals[j]*sol1[aj[j]];
//      }
//    }
//    for(j=ai[i]; j<nz; j++) {
//      b02[i]+=vals[j]*sol1[aj[j]];
//    }
    spec48_rpesol_(insize+j1*insizes,irnereg[j2],valereg[j2],keepreg[j2],b02,b03,cntl,rinfo,error1,icntl,info,w,iw);
    free(b02);
    for(j=0; j<nrow; j++)sol2[j]-=b03[j];
    //free(b03);
    sol2+=nrow;
    //b03+=nrow;
  }
  //printf("OKKK111222\n");
  return true;
}

bool NDBBD_sol_nread1(PetscInt rank, int begmat,int nreg,int * insize,int insizes, Mat **submatCij,Mat **submatBij,ha_cgetype *b,ha_cgetype *sol,int** irnereg,int** keepreg,ha_cgetype** valereg,ha_cgetype *cntl,ha_cgetype *rinfo,ha_cgetype *error1,int *icntl,int *info,ha_cgetype *w,int *iw,ha_cgetype *b02) {
  //FILE* fp1,*fp2,*fp3;
  PetscScalar *vals;//,*valsc;//,vecval;
  PetscInt *ai,*aj,nrow,nz;//,ncol,*aic,*ajc
  PetscInt i,j,j1,j2;//,indx01,la1;
  ha_cgetype *b01,*b03,*sol1,*sol2;//,*b02
  int * insize1;
  b01=b;
  sol1=sol;
  //printf("OKKK111222\n");
  for (j1=begmat; j1<nreg+begmat; j1++) {
    j2=j1-begmat;//-begmat;
    j=j1*insizes;
    //printf("j2 %d\n",j2);
    insize1=insize+j;
    insize1[16]=ceil((insize1[9]/100.0)*insize1[2]);
    spec48m_rpesol_(insize1,irnereg[j2],valereg[j2],keepreg[j2],b01,sol1,cntl,rinfo,error1,icntl,info,w,iw);//insize+j
    b01+=insize[j];
    sol1+=insize[j];
  }
  //printf("OKKK111222\n");
  sol2=sol;
  for (j1=begmat; j1<nreg+begmat; j1++) {
    Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatBij[j1][0]->data;
    ai= ab->i;
    aj= ab->j;
    vals=ab->a;
    nz=ab->nz;
    nrow=submatBij[j1][0]->rmap->n;
    spar_mulmin_(sol2,&nrow,&nz,ai,aj,vals,b01);
//    for(i=0; i<nrow-1; i++) {
//      for(j=ai[i]; j<ai[i+1]; j++) {
//        b01[i]-=vals[j]*sol2[aj[j]];
//      }
//    }
//    for(j=ai[i]; j<nz; j++) {
//      b01[i]-=vals[j]*sol2[aj[j]];
//    }
    sol2+=insize[j1*insizes];
  }
  j2=j1-begmat;//-begmat;
  insize1=insize+j1*insizes;
  insize1[16]=ceil((insize1[9]/100.0)*insize1[2]);
  spec48m_rpesol_(insize1,irnereg[j2],valereg[j2],keepreg[j2],b01,sol1,cntl,rinfo,error1,icntl,info,w,iw);//insize+j1*insizes
  sol2=sol;
  //printf("OKKK111222\n");
  b03=b;
  for (j1=begmat; j1<nreg+begmat; j1++) {
    j2=j1-begmat;//-begmat;
    Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatCij[j1][0]->data;//*aa=subA->data;
    ai= ac->i;
    aj= ac->j;
    vals=ac->a;
    nz=ac->nz;
    nrow=submatCij[j1][0]->rmap->n;
    //ha_cgetype *b02 = (ha_cgetype*)calloc(nrow,sizeof(ha_cgetype));
    //ha_cgetype *b03 = (ha_cgetype*)calloc(nrow,sizeof(ha_cgetype));
    spar_mulnoadd_(sol1,&nrow,&nz,ai,aj,vals,b02);
//    for(i=0; i<nrow-1; i++) {
//      for(j=ai[i]; j<ai[i+1]; j++) {
//        b02[i]+=vals[j]*sol1[aj[j]];
//      }
//    }
//    for(j=ai[i]; j<nz; j++) {
//      b02[i]+=vals[j]*sol1[aj[j]];
//    }
    insize1=insize+j1*insizes;
    insize1[16]=ceil((insize1[9]/100.0)*insize1[2]);
    spec48m_rpesol_(insize1,irnereg[j2],valereg[j2],keepreg[j2],b02,b03,cntl,rinfo,error1,icntl,info,w,iw);//insize+j1*insizes
    //free(b02);
    for(j=0; j<nrow; j++)sol2[j]-=b03[j];
    //free(b03);
    sol2+=nrow;
    //b03+=nrow;
  }
  //printf("OKKK111222\n");
  return true;
}

int HaReduce(ha_cgetype *vecbivi,forint vecbivisize,PetscInt mpisize,PetscInt rank,PetscInt targetrank) {
  if(mpisize==1)return 0;
  MPI_Status   status;
  int i,j,j1,j2,j3;
  //ha_cgetype *vecbivi0;
  j=0;
  if(rank!=targetrank){
  #pragma omp parallel private(i) reduction(+:j)
  {
  j=0;
  #pragma omp for
    for(i=0; i<vecbivisize; i++)if(vecbivi[i]!=0)j++;
  }
  }
  if(rank==targetrank) j=0;//(VecSize-sumrowcolin)*(VecSize-sumrowcolin);
  //printf("jjjju %d rank %d vs %d\n",j,rank,vecbivisize);
  int *biviindx= (int *) calloc (j,sizeof(int));
  ha_cgetype *vecbivi0= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
  j2=0;
  //printf("jjjju1 %d rank %d\n",j,rank);
  if(rank!=targetrank)for(i=0; i<vecbivisize; i++) {
      if(vecbivi[i]!=0) {
        vecbivi[j2]=vecbivi[i];
        biviindx[j2]=i;
        j2++;
      }
      //printf("i %d\n",i);
    }
  //MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Reduce(&j,&j2,1, MPI_INT, MPI_MAX,targetrank,PETSC_COMM_WORLD);
  //printf("j2maxu %d\n",j2);
  if(rank==targetrank) {
    vecbivi0=realloc(vecbivi0,j2*sizeof(ha_cgetype));
    biviindx=realloc(biviindx,j2*sizeof(int));
  }
  //if(rank!=0)vecbivi=realloc(vecbivi,j*sizeof(ha_cgetype));
  //printf("jjjju1 %d rank %d\n",j,rank);

  for(j1=0; j1<targetrank; j1++) {
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,j1, PETSC_COMM_WORLD);
    j3=j1;
    MPI_Bcast(&j3,1, MPI_INT,j1, PETSC_COMM_WORLD);
    if(rank==j1) {
      if(SORD==1)MPI_Send(vecbivi,j, MPI_DOUBLE, targetrank, 10, PETSC_COMM_WORLD);
      else MPI_Send(vecbivi,j, MPI_FLOAT, targetrank, 10, PETSC_COMM_WORLD);
      MPI_Send(biviindx,j, MPI_INT, targetrank, 11, PETSC_COMM_WORLD);
    }
    //printf("j2 %d\n",j2);
    if(rank==targetrank) {
      //vecbivi0=realloc(vecbivi0,j2*sizeof(ha_cgetype));
      //biviindx=realloc(biviindx,j2*sizeof(int));
      //printf("j3 %d\n",j3);
      if(SORD==1)MPI_Recv(vecbivi0,j2, MPI_DOUBLE, j3, 10, PETSC_COMM_WORLD,&status);
      else MPI_Recv(vecbivi0,j2, MPI_FLOAT, j3, 10, PETSC_COMM_WORLD,&status);
      MPI_Recv(biviindx,j2, MPI_INT, j3, 11, PETSC_COMM_WORLD,&status);
        #pragma omp parallel private(i)
        {
        #pragma omp for
      for(i=0; i<j2; i++) {
        //printf("biviindx %d\n",biviindx[i]);
        #pragma omp atomic
        vecbivi[biviindx[i]]+=vecbivi0[i];
      }
        }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  for(j1=targetrank+1; j1<mpisize; j1++) {
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,j1, PETSC_COMM_WORLD);
    j3=j1;
    MPI_Bcast(&j3,1, MPI_INT,j1, PETSC_COMM_WORLD);
    if(rank==j1) {
      if(SORD==1)MPI_Send(vecbivi,j, MPI_DOUBLE, targetrank, 10, PETSC_COMM_WORLD);
      else MPI_Send(vecbivi,j, MPI_FLOAT, targetrank, 10, PETSC_COMM_WORLD);
      MPI_Send(biviindx,j, MPI_INT, targetrank, 11, PETSC_COMM_WORLD);
    }
    //printf("j2 %d\n",j2);
    if(rank==targetrank) {
      //vecbivi0=realloc(vecbivi0,j2*sizeof(ha_cgetype));
      //biviindx=realloc(biviindx,j2*sizeof(int));
      //printf("j3 %d\n",j3);
      if(SORD==1)MPI_Recv(vecbivi0,j2, MPI_DOUBLE, j3, 10, PETSC_COMM_WORLD,&status);
      else MPI_Recv(vecbivi0,j2, MPI_FLOAT, j3, 10, PETSC_COMM_WORLD,&status);
      MPI_Recv(biviindx,j2, MPI_INT, j3, 11, PETSC_COMM_WORLD,&status);
        #pragma omp parallel private(i)
        {
        #pragma omp for
      for(i=0; i<j2; i++) {
        //printf("biviindx %d\n",biviindx[i]);
        #pragma omp atomic
        vecbivi[biviindx[i]]+=vecbivi0[i];
      }
        }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  free(biviindx);
  free(vecbivi0);
  return 1;
}

int HaReduceNoComp(ha_cgetype *vecbivi,forint vecbivisize,PetscInt mpisize,PetscInt rank,PetscInt targetrank) {
  if(mpisize==1)return 0;
  MPI_Status   status;
  int i,j,j1,j2,j3;
  ha_cgetype *vecbivi0= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));

  if(rank==targetrank) {
    vecbivi0=realloc(vecbivi0,vecbivisize*sizeof(ha_cgetype));
  }

  for(j1=0; j1<targetrank; j1++) {
    if(rank==j1) {
      if(SORD==1)MPI_Send(vecbivi,vecbivisize, MPI_DOUBLE, targetrank, 101, PETSC_COMM_WORLD);
      else MPI_Send(vecbivi,vecbivisize, MPI_FLOAT, targetrank, 101, PETSC_COMM_WORLD);
    }
    if(rank==targetrank) {
      if(SORD==1)MPI_Recv(vecbivi0,vecbivisize, MPI_DOUBLE, j1, 101, PETSC_COMM_WORLD,&status);
      else MPI_Recv(vecbivi0,vecbivisize, MPI_FLOAT, j1, 101, PETSC_COMM_WORLD,&status);
        #pragma omp parallel private(i)
        {
        #pragma omp for
      for(i=0; i<vecbivisize; i++) {
        #pragma omp atomic
        vecbivi[i]+=vecbivi0[i];
      }
        }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  for(j1=targetrank+1; j1<mpisize; j1++) {
    if(rank==j1) {
      if(SORD==1)MPI_Send(vecbivi,vecbivisize, MPI_DOUBLE, targetrank, 102, PETSC_COMM_WORLD);
      else MPI_Send(vecbivi,vecbivisize, MPI_FLOAT, targetrank, 102, PETSC_COMM_WORLD);
    }
    if(rank==targetrank) {
      if(SORD==1)MPI_Recv(vecbivi0,vecbivisize, MPI_DOUBLE, j1, 102, PETSC_COMM_WORLD,&status);
      else MPI_Recv(vecbivi0,vecbivisize, MPI_FLOAT, j1, 102, PETSC_COMM_WORLD,&status);
        #pragma omp parallel private(i)
        {
        #pragma omp for
      for(i=0; i<vecbivisize; i++) {
        #pragma omp atomic
        vecbivi[i]+=vecbivi0[i];
      }
        }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  free(vecbivi0);
  return 1;
}

/*int HaReduceNoComp(ha_cgetype *vecbivi,forint vecbivisize,PetscInt mpisize,PetscInt rank,PetscInt targetrank) {
  if(mpisize==1)return 0;
  MPI_Status   status;
  int i,j,j1,j2,j3;
  forint halfvec=vecbivisize/2;
  ha_cgetype *vecbivi0= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));

  if(rank==targetrank) {
    vecbivi0=realloc(vecbivi0,(vecbivisize-halfvec)*sizeof(ha_cgetype));
  }

  for(j1=0; j1<targetrank; j1++) {
    if(rank==j1) {
      MPI_Send(vecbivi+halfvec,vecbivisize-halfvec, MPI_DOUBLE, targetrank, 10, PETSC_COMM_WORLD);
      vecbivi=realloc(vecbivi,halfvec*sizeof(ha_cgetype));
      vecbivi0=realloc(vecbivi0,halfvec*sizeof(ha_cgetype));
      MPI_Recv(vecbivi0,halfvec, MPI_DOUBLE, targetrank, 15, PETSC_COMM_WORLD,&status);
    }
    if(rank==targetrank) {
      MPI_Recv(vecbivi0,vecbivisize-halfvec, MPI_DOUBLE, j1, 10, PETSC_COMM_WORLD,&status);
      MPI_Send(vecbivi,halfvec, MPI_DOUBLE, j1, 15, PETSC_COMM_WORLD);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    if(rank==j1) {
      for(i=0; i<halfvec; i++) {
        vecbivi[i]+=vecbivi0[i];
      }
      MPI_Send(vecbivi,halfvec, MPI_DOUBLE, targetrank, 17, PETSC_COMM_WORLD);
      vecbivi=realloc(vecbivi,sizeof(ha_cgetype));
      vecbivi0=realloc(vecbivi0,sizeof(ha_cgetype));
    }
    if(rank==targetrank) {
      for(i=0; i<vecbivisize-halfvec; i++) {
        vecbivi[i+halfvec]+=vecbivi0[i];
      }
      MPI_Recv(vecbivi0,halfvec, MPI_DOUBLE, j1, 17, PETSC_COMM_WORLD,&status);
      for(i=0; i<halfvec; i++) {
        vecbivi[i]+=vecbivi0[i];
      }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }

  for(j1=targetrank+1; j1<mpisize; j1++) {
    if(rank==j1) {
      MPI_Send(vecbivi+halfvec,vecbivisize-halfvec, MPI_DOUBLE, targetrank, 11, PETSC_COMM_WORLD);
      vecbivi=realloc(vecbivi,halfvec*sizeof(ha_cgetype));
      vecbivi0=realloc(vecbivi0,halfvec*sizeof(ha_cgetype));
      MPI_Recv(vecbivi0,halfvec, MPI_DOUBLE, targetrank, 16, PETSC_COMM_WORLD,&status);
    }
    if(rank==targetrank) {
      MPI_Recv(vecbivi0,vecbivisize-halfvec, MPI_DOUBLE, j1, 11, PETSC_COMM_WORLD,&status);
      MPI_Send(vecbivi,halfvec, MPI_DOUBLE, j1, 16, PETSC_COMM_WORLD);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    if(rank==j1) {
      for(i=0; i<halfvec; i++) {
        vecbivi[i]+=vecbivi0[i];
      }
      MPI_Send(vecbivi,halfvec, MPI_DOUBLE, targetrank, 18, PETSC_COMM_WORLD);
      vecbivi=realloc(vecbivi,sizeof(ha_cgetype));
      vecbivi0=realloc(vecbivi0,sizeof(ha_cgetype));
    }
    if(rank==targetrank) {
      for(i=0; i<vecbivisize-halfvec; i++) {
        vecbivi[i+halfvec]+=vecbivi0[i];
      }
      MPI_Recv(vecbivi0,halfvec, MPI_DOUBLE, j1, 18, PETSC_COMM_WORLD,&status);
      for(i=0; i<halfvec; i++) {
        vecbivi[i]+=vecbivi0[i];
      }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  free(vecbivi0);
  return 1;
}*/

int HaReducebivi(ha_cgetype *vecbivi,int *nz,int *biviindx1,int *bivinzrow,int *bivinzcol,int *bivinzrow0,int *bivinzcol0,int *bivinzrow1,int *bivinzcol1,int bivirowsize,int bivicolsize,PetscInt Dcolsize,PetscInt mpisize,PetscInt rank,PetscInt targetrank) {
  if(mpisize==1)return 0;
  MPI_Status   status;
  int i,j,j1,j2,j3,j4,vecbivisize=bivirowsize*bivicolsize,nz0,nz1,nz2,gap;
  //int bivicolsize0=bivicolsize;
  //printf("jjjj1 %d rank %d\n",j,rank);
  //MPI_Bcast(&bivicolsize0,sizeof(int), MPI_BYTE,targetrank, PETSC_COMM_WORLD);
  //printf("jjjj1 %d rank %d\n",j,rank);
  //ha_cgetype *vecbivi0;
  j=0;
  //if(rank!=targetrank)
  for(i=0; i<vecbivisize; i++)if(vecbivi[i]!=0)j++;
  //if(rank==targetrank) j=0;//(VecSize-sumrowcolin)*(VecSize-sumrowcolin);
  //printf("jjjj %d rank %d\n",j,rank);
  //int *biviindx= (int *) calloc (j,sizeof(int));
  biviindx1=realloc(biviindx1,j*sizeof(int));
  int *biviindx0= (int *) calloc (1,sizeof(int));
  ha_cgetype *vecbivi0= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
  j2=0;
  //printf("jjjj1 %d rank %d col %d col0 %d\n",j,rank,bivicolsize,bivicolsize0);
  //if(rank!=targetrank)
  {
    for(i=0; i<vecbivisize; i++) {
      if(vecbivi[i]!=0) {
        vecbivi[j2]=vecbivi[i];
        j1=i/bivicolsize;
        j3=i-j1*bivicolsize;
        //printf("j1 %d j3 %d\n",i,bivinzcol[bivinzcol1[j3]]-1+(bivinzrow[bivinzrow1[j1]]-1)*bivicolsize);
        //biviindx[j2]=bivinzcol0[bivinzcol1[j3]]-1+(bivinzrow0[bivinzrow1[j1]]-1)*bivicolsize0;
        biviindx1[j2]=bivinzcol1[j3]+bivinzrow1[j1]*Dcolsize;
        //if(rank==0)printf("v %d\n",biviindx[j2]);
        j2++;
      }
    }
    vecbivi=realloc(vecbivi,j2*sizeof(ha_cgetype));
  }
  nz[0]=j2;
  printf("jjjj1 %d j2 %d rank %d\n",j,j2,rank);
  MPI_Reduce(&j,&j2,1, MPI_INT, MPI_MAX,targetrank,PETSC_COMM_WORLD);
  if(rank==targetrank) {
    vecbivi0=realloc(vecbivi0,j2*sizeof(ha_cgetype));
    biviindx0=realloc(biviindx0,j2*sizeof(int));
    printf("j2max %d\n",j2);
  }
  //if(rank!=0)vecbivi=realloc(vecbivi,j*sizeof(ha_cgetype));

  for(j1=0; j1<targetrank; j1++) {
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,j1, PETSC_COMM_WORLD);
    j3=j1;
    MPI_Bcast(&j3,1, MPI_INT,j1, PETSC_COMM_WORLD);
    if(rank==j1) {
      if(SORD==1)MPI_Send(vecbivi,j, MPI_DOUBLE, targetrank, 10, PETSC_COMM_WORLD);
      else MPI_Send(vecbivi,j, MPI_FLOAT, targetrank, 10, PETSC_COMM_WORLD);
      MPI_Send(biviindx1,j, MPI_INT, targetrank, 11, PETSC_COMM_WORLD);
    }
    printf("j2 %d\n",j2);
    if(rank==targetrank) {
      //vecbivi0=realloc(vecbivi0,j2*sizeof(ha_cgetype));
      //biviindx=realloc(biviindx,j2*sizeof(int));
      //printf("j3 %d\n",j3);
      if(SORD==1)MPI_Recv(vecbivi0,j2, MPI_DOUBLE, j3, 10, PETSC_COMM_WORLD,&status);
      else MPI_Recv(vecbivi0,j2, MPI_FLOAT, j3, 10, PETSC_COMM_WORLD,&status);
      MPI_Recv(biviindx0,j2, MPI_INT, j3, 11, PETSC_COMM_WORLD,&status);
      printf("in %d in0 %d j1 %d\n",biviindx1[0],biviindx0[0],j1);


      j=0;
      nz1=0;
      nz2=nz[0];
      printf("nz2 %d\n",nz2);
      j4=j2-1;
      printf("in %d in0 %d\n",biviindx1[0],biviindx0[0]);
      for(i=0; i<nz2; i++) {
        //printf("j %d i %d\n",j,i);
        if(biviindx1[i]<biviindx0[j]) {
          nz1++;
          continue;
        }
        if(biviindx1[i]==biviindx0[j]) {
          nz1++;
          j++;
          if(j==j4)break;
          continue;
        }
        for(j3=j; j3<j2; j3++) {
          //printf("j3 %d\n",j3);
          if(biviindx1[i]<biviindx0[j3]) {
            nz1++;
            break;
          }
          if(biviindx1[i]==biviindx0[j3]) {
            nz1++;
            j3++;
            break;
          }
          nz1++;
        }
        j=j3;
        //nz1++;
        if(j==j4)break;
      }
      printf("i %d j %d nz1 %d\n",i,j,nz1);
      if(i<nz2-1)nz1+=nz2-1-i;
      if(j<j4)nz1+=j4-j;
      printf("nz1 %d\n",nz1);
      vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
      //printf("hare nz1 %d\n",nz1);
      //if(vecbivi==NULL)printf("Relloc failed!");
      biviindx1=realloc(biviindx1,nz1*sizeof(int));
      //if(biviindx==NULL)printf("Relloc failed!");
      //printf("hare nz1 %d\n",nz1);
      j=j2-1;
      j1=nz1-1;
      for(i=nz2-1; i>-1; i--) {
        //printf("i %d j1 %d j %d v %d v0 %d v1 %d\n",i,j1,j,biviindx[i],biviindx0[j],biviindx[j1+1]);
        if(biviindx1[i]>biviindx0[j]) {
          vecbivi[j1]=vecbivi[i];
          biviindx1[j1]=biviindx1[i];
          j1--;
          continue;
        }
        if(biviindx1[i]==biviindx0[j]) {
          vecbivi[j1]=vecbivi[i]+vecbivi0[j];
          biviindx1[j1]=biviindx1[i];
          j--;
          j1--;
          if(j==-1)break;
          continue;
        }
        for(j3=j; j3>-1; j3--) {
          if(biviindx1[i]>biviindx0[j3]) {
            vecbivi[j1]=vecbivi[i];
            biviindx1[j1]=biviindx1[i];
            j1--;
            break;
          }
          if(biviindx1[i]==biviindx0[j3]) {
            vecbivi[j1]=vecbivi[i]+vecbivi0[j3];
            biviindx1[j1]=biviindx1[i];
            j3--;
            j1--;
            break;
          }
          biviindx1[j1]=biviindx0[j3];
          vecbivi[j1]=vecbivi0[j3];
          j1--;
        }
        j=j3;
        //biviindx[j1]=biviindx[i];
        //vecbivi[j1]=vecbivi[i];
        //j1--;
        if(j==-1)break;
      }
      printf("j1 %d j %d i %d v %d v0 %d\n",j1,j,i,biviindx1[i+1],biviindx0[j+1]);
      if(j>-1)for(i=j; i>-1; i--) {
          biviindx1[j1]=biviindx0[i];
          vecbivi[j1]=vecbivi0[i];
          j1--;
        }
      printf("j1 %d j %d i %d vm1 %d\n",j1,j,i,biviindx1[j1+1]);
      for(i=0; i<2; i++)printf("vout %d\n",biviindx1[i]);
      nz[0]=nz1;


      //HaSpaVecAdd(vecbivi,biviindx1,nz,vecbivi0,biviindx0,j2);
      printf("out in %d in0 %d j1 %d\n",biviindx1[0],biviindx0[0],j1);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  printf("here!!!\n");
  for(j1=targetrank+1; j1<mpisize; j1++) {
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,j1, PETSC_COMM_WORLD);
    j3=j1;
    MPI_Bcast(&j3,1, MPI_INT,j1, PETSC_COMM_WORLD);
    if(rank==j1) {
      if(SORD==1)MPI_Send(vecbivi,j, MPI_DOUBLE, targetrank, 10, PETSC_COMM_WORLD);
      else MPI_Send(vecbivi,j, MPI_FLOAT, targetrank, 10, PETSC_COMM_WORLD);
      MPI_Send(biviindx1,j, MPI_INT, targetrank, 11, PETSC_COMM_WORLD);
    }
    //printf("j2 %d\n",j2);
    if(rank==targetrank) {
      //vecbivi0=realloc(vecbivi0,j2*sizeof(ha_cgetype));
      //biviindx=realloc(biviindx,j2*sizeof(int));
      //printf("j3 %d\n",j3);
      if(SORD==1)MPI_Recv(vecbivi0,j2, MPI_DOUBLE, j3, 10, PETSC_COMM_WORLD,&status);
      else MPI_Recv(vecbivi0,j2, MPI_FLOAT, j3, 10, PETSC_COMM_WORLD,&status);
      MPI_Recv(biviindx0,j2, MPI_INT, j3, 11, PETSC_COMM_WORLD,&status);
      HaSpaVecAdd(vecbivi,biviindx1,nz,vecbivi0,biviindx0,j2);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  free(biviindx0);
  free(vecbivi0);
  printf("here1!!!\n");
  printf("out in1 %d\n",biviindx1[0]);
  //biviindx1=biviindx;
  //nz[0]=nz0;
  return 1;
}

int HaSpaVecAdd(ha_cgetype *vecbivi,int *biviindx,int *nz0,ha_cgetype *vecbivi0,int *biviindx0,int j2) {
  int j,j1,j3,j4,nz1,nz2,i;
  j=0;
  nz1=0;
  nz2=nz0[0];
  printf("nz2 %d\n",nz2);
  j4=j2-1;
  printf("in %d in0 %d\n",biviindx[0],biviindx0[0]);
  for(i=0; i<nz2; i++) {
    //printf("j %d i %d\n",j,i);
    if(biviindx[i]<biviindx0[j]) {
      nz1++;
      continue;
    }
    if(biviindx[i]==biviindx0[j]) {
      nz1++;
      j++;
      if(j==j4)break;
      continue;
    }
    for(j3=j; j3<j2; j3++) {
      //printf("j3 %d\n",j3);
      if(biviindx[i]<biviindx0[j3]) {
        nz1++;
        break;
      }
      if(biviindx[i]==biviindx0[j3]) {
        nz1++;
        j3++;
        break;
      }
      nz1++;
    }
    j=j3;
    //nz1++;
    if(j==j4)break;
  }
  printf("i %d j %d nz1 %d\n",i,j,nz1);
  if(i<nz2-1)nz1+=nz2-1-i;
  if(j<j4)nz1+=j4-j;
  printf("nz1 %d\n",nz1);
  vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
  //printf("hare nz1 %d\n",nz1);
  //if(vecbivi==NULL)printf("Relloc failed!");
  biviindx=realloc(biviindx,nz1*sizeof(int));
  //if(biviindx==NULL)printf("Relloc failed!");
  //printf("hare nz1 %d\n",nz1);
  j=j2-1;
  j1=nz1-1;
  for(i=nz2-1; i>-1; i--) {
    //printf("i %d j1 %d j %d v %d v0 %d v1 %d\n",i,j1,j,biviindx[i],biviindx0[j],biviindx[j1+1]);
    if(biviindx[i]>biviindx0[j]) {
      vecbivi[j1]=vecbivi[i];
      biviindx[j1]=biviindx[i];
      j1--;
      continue;
    }
    if(biviindx[i]==biviindx0[j]) {
      vecbivi[j1]=vecbivi[i]+vecbivi0[j];
      biviindx[j1]=biviindx[i];
      j--;
      j1--;
      if(j==-1)break;
      continue;
    }
    for(j3=j; j3>-1; j3--) {
      if(biviindx[i]>biviindx0[j3]) {
        vecbivi[j1]=vecbivi[i];
        biviindx[j1]=biviindx[i];
        j1--;
        break;
      }
      if(biviindx[i]==biviindx0[j3]) {
        vecbivi[j1]=vecbivi[i]+vecbivi0[j3];
        biviindx[j1]=biviindx[i];
        j3--;
        j1--;
        break;
      }
      biviindx[j1]=biviindx0[j3];
      vecbivi[j1]=vecbivi0[j3];
      j1--;
    }
    j=j3;
    //biviindx[j1]=biviindx[i];
    //vecbivi[j1]=vecbivi[i];
    //j1--;
    if(j==-1)break;
  }
  printf("j1 %d j %d i %d v %d v0 %d\n",j1,j,i,biviindx[i+1],biviindx0[j+1]);
  if(j>-1)for(i=j; i>-1; i--) {
      biviindx[j1]=biviindx0[i];
      vecbivi[j1]=vecbivi0[i];
      j1--;
    }
  printf("j1 %d j %d i %d vm1 %d\n",j1,j,i,biviindx[j1+1]);
  for(i=0; i<2; i++)printf("vout %d\n",biviindx[i]);
  nz0[0]=nz1;
}
/*int HaMatGetSubMats(Mat A,uvadd indxbsize,PetscInt *indicesB,uvadd *ha_ndblocks,uvadd *begblock, uvadd *countvarintra1, uvadd *ha_cols,PetscInt mpisize,PetscInt rank,forint nmatin,PetscInt origrank,int **nzbs,int **bns,int **bms,int **irnbs,int **jcnbs,PetscScalar **valbs) {
  int bfirst,bend,i,j,j0,j1,j2,j3,j4;
  int nz,nz1;
  int j5;
  PetscInt *ai,*aj;
  PetscScalar *vals;
  Mat *submatB,submatBT;
  IS *colindices,*rowindices;
  MPI_Status   status;
  i=0;
  for(i=0; i<mpisize-1; i++) {
    j0=nmatin;
    MPI_Bcast(&j0,sizeof(j0), MPI_BYTE,i, PETSC_COMM_WORLD);
    //printf("i0 %d rank %d\n",i,rank);
    if(rank==mpisize-1) {
      //IS *colindices= (IS *) calloc (j0,sizeof(**IS));
      PetscMalloc(j0*sizeof(IS **),&colindices);
      PetscMalloc(j0*sizeof(IS **),&rowindices);
      for(j=0; j<j0; j++) {
        PetscInt *indices= (PetscInt *) calloc (ha_ndblocks[j+begblock[i]],sizeof(PetscInt));
        bfirst=countvarintra1[j+begblock[i]];
        for(j1=0; j1<ha_ndblocks[j+begblock[i]]; j1++) {
          indices[j1]=ha_cols[bfirst+j1]+bfirst;
          //printf("indice %d j1 %d size %d\n",indices[j1],j1,ha_ndblocks[j+begblock[i]]);
        }
        ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[j+begblock[i]],indices,PETSC_COPY_VALUES,colindices+j);
        //ISView(rowindices[j],PETSC_VIEWER_STDOUT_SELF);
        free(indices);
        ISCreateGeneral(PETSC_COMM_SELF,indxbsize,indicesB,PETSC_COPY_VALUES,rowindices+j);
      }
      PetscMalloc(j0*sizeof(Mat *),&submatB);
      //printf("j0 %d\n",j0);
    }
    j1=j0;
    if(rank!=mpisize-1)j1=0;
    MatCreateSubMatrices(A,j1,rowindices,colindices,MAT_INITIAL_MATRIX,&submatB);
    if(rank==mpisize-1) {
      for(j=0; j<j0; j++) {
        ISDestroy(colindices+j);
        ISDestroy(rowindices+j);
      }
      PetscFree(colindices);
      PetscFree(rowindices);
      //MatView(submatB[0],0);
      nz=0;
      nz1=0;
      for(j=0; j<j0; j++) {
        Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatB[j]->data;
        //nzbs[j]=ab->nz;
        nz+=ab->nz;//nzbs[j];
        nz1+=submatB[j]->cmap->n;//nzbs[j];
      }
      if(i<mpisize-1) {
        MPI_Send(&nz,1, MPI_INT, i, 0, PETSC_COMM_WORLD);
        MPI_Send(&nz1,1, MPI_INT, i, 7, PETSC_COMM_WORLD);
      }
      //printf("i %d rank %d nz %d\n",i,rank,nz);
    }
    if(rank==i) {
      MPI_Recv(&nz,1, MPI_INT, mpisize-1, 0, PETSC_COMM_WORLD,&status);
      MPI_Recv(&nz1,1, MPI_INT, mpisize-1, 7, PETSC_COMM_WORLD,&status);
    }
    if(rank==i||rank==mpisize-1) {
      //*nzbs= calloc (*nzbs,j0*sizeof(int));
      //*bns= calloc (*bns,j0*sizeof(int));
      //*bms= calloc (*bms,j0*sizeof(int));
      //*jcnbs= calloc (*jcnbs,nz*sizeof(int));
      //*irnbs= calloc (*irnbs,nz1*sizeof(int));
      //*valbs= calloc (*valbs,nz*sizeof(PetscScalar));
      *nzbs= (int*)calloc (j0,sizeof(int));
      *bns= (int*)calloc (j0,sizeof(int));
      *bms= (int*)calloc (j0,sizeof(int));
      *jcnbs= (int*)calloc (nz,sizeof(int));
      *irnbs= (int*)calloc (nz1,sizeof(int));
      *valbs= (PetscScalar*)calloc (nz,sizeof(PetscScalar));
      //printf("i5 %d rank %d nz %d\n",i,rank,nz);
    }
    if(rank==mpisize-1) {
      j3=0;
      j4=0;
      for(j=0; j<j0; j++) {
        MatTranspose(submatB[j],MAT_INITIAL_MATRIX,&submatBT);
        MatDestroy_SeqAIJ(submatB[j]);
        Mat_SeqAIJ         *abt=(Mat_SeqAIJ*)submatBT->data;
        ai= abt->i;
        aj= abt->j;
        vals=abt->a;
        (*nzbs)[j]=abt->nz;
        (*bns)[j]=submatBT->rmap->n;
        (*bms)[j]=submatBT->cmap->n;
        //printf("j %d btn %d bms %d nz %d\n",j,submatBT->cmap->n,(*bms)[j],(*nzbs)[j]);
        for(j1=0; j1<(*bns)[j]; j1++) {
          (*irnbs)[j4]=ai[j1];
          j4++;
        }
        for(j1=0; j1<(*nzbs)[j]; j1++) {
          (*jcnbs)[j3]=aj[j1];
          //printf("jcn %d\n",(*jcnbs)[j3]);
          (*valbs)[j3]=vals[j1];
          j3++;
        }
        MatDestroy_SeqAIJ(submatBT);
      }
      free(submatB);
    }
    if(i<mpisize-1) {
      if(rank==mpisize-1) {
        MPI_Send(*valbs,nz, MPI_DOUBLE,i, 31, PETSC_COMM_WORLD);
        MPI_Send(*irnbs,nz1*sizeof(int), MPI_BYTE,i, 32, PETSC_COMM_WORLD);
        MPI_Send(*jcnbs,nz*sizeof(int), MPI_BYTE,i, 33, PETSC_COMM_WORLD);
        MPI_Send(*bns,j0, MPI_INT,i, 34, PETSC_COMM_WORLD);
        MPI_Send(*bms,j0, MPI_INT,i, 35, PETSC_COMM_WORLD);
        MPI_Send(*nzbs,j0*sizeof(int), MPI_BYTE,i, 36, PETSC_COMM_WORLD);
      }
      if(rank==i) {
        MPI_Recv(*valbs,nz, MPI_DOUBLE,mpisize-1, 31, PETSC_COMM_WORLD,&status);
        MPI_Recv(*irnbs,nz1*sizeof(int), MPI_BYTE,mpisize-1, 32, PETSC_COMM_WORLD,&status);
        MPI_Recv(*jcnbs,nz*sizeof(int), MPI_BYTE,mpisize-1, 33, PETSC_COMM_WORLD,&status);
        MPI_Recv(*bns,j0, MPI_INT,mpisize-1, 34, PETSC_COMM_WORLD,&status);
        MPI_Recv(*bms,j0, MPI_INT,mpisize-1, 35, PETSC_COMM_WORLD,&status);
        MPI_Recv(*nzbs,j0*sizeof(int), MPI_BYTE,mpisize-1, 36, PETSC_COMM_WORLD,&status);
        //for(j1=0; j1<j0; j1++)printf("nzzzzzz %d rank %d\n",(*nzbs)[j1],rank);
      }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    if(rank==0)for(j1=0; j1<j0; j1++)printf("nzzzzzz %d rank %d i %d\n",(*nzbs)[j1],rank,i);
  }

  i=mpisize-1;
  j0=nmatin;
  if(rank==mpisize-1) {
    //IS *colindices= (IS *) calloc (j0,sizeof(**IS));
    PetscMalloc(j0*sizeof(IS **),&colindices);
    PetscMalloc(j0*sizeof(IS **),&rowindices);
    for(j=0; j<j0; j++) {
      PetscInt *indices= (PetscInt *) calloc (ha_ndblocks[j+begblock[i]],sizeof(PetscInt));
      bfirst=countvarintra1[j+begblock[i]];
      for(j1=0; j1<ha_ndblocks[j+begblock[i]]; j1++) {
        indices[j1]=ha_cols[bfirst+j1]+bfirst;
        //printf("indice %d j1 %d size %d\n",indices[j1],j1,ha_ndblocks[j+begblock[i]]);
      }
      ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[j+begblock[i]],indices,PETSC_COPY_VALUES,colindices+j);
      //ISView(rowindices[j],PETSC_VIEWER_STDOUT_SELF);
      free(indices);
      ISCreateGeneral(PETSC_COMM_SELF,indxbsize,indicesB,PETSC_COPY_VALUES,rowindices+j);
    }
    PetscMalloc(j0*sizeof(Mat *),&submatB);
    //printf("j0 %d\n",j0);
  }
  j1=j0;
  if(rank!=mpisize-1)j1=0;
  MatCreateSubMatrices(A,j1,rowindices,colindices,MAT_INITIAL_MATRIX,&submatB);
  if(rank==mpisize-1) {
    for(j=0; j<j0; j++) {
      ISDestroy(colindices+j);
      ISDestroy(rowindices+j);
    }
    PetscFree(colindices);
    PetscFree(rowindices);
    nz=0;
    nz1=0;
    for(j=0; j<j0; j++) {
      Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatB[j]->data;
      //nzbs[j]=ab->nz;
      nz+=ab->nz;//nzbs[j];
      nz1+=submatB[j]->cmap->n;//nzbs[j];
    }
    //*nzbs= realloc (*nzbs,j0*sizeof(int));
    //*bns= realloc (*bns,j0*sizeof(int));
    //*bms= realloc (*bms,j0*sizeof(int));
    //*jcnbs= realloc (*jcnbs,nz*sizeof(int));
    //*irnbs= realloc (*irnbs,nz1*sizeof(int));
    //*valbs= realloc (*valbs,nz*sizeof(PetscScalar));
    *nzbs= (int*)calloc (j0,sizeof(int));
    *bns= (int*)calloc (j0,sizeof(int));
    *bms= (int*)calloc (j0,sizeof(int));
    *jcnbs= (int*)calloc (nz,sizeof(int));
    *irnbs= (int*)calloc (nz1,sizeof(int));
    *valbs= (PetscScalar*)calloc (nz,sizeof(PetscScalar));
    j3=0;
    j4=0;
    for(j=0; j<j0; j++) {
      MatTranspose(submatB[j],MAT_INITIAL_MATRIX,&submatBT);
      MatDestroy_SeqAIJ(submatB[j]);
      Mat_SeqAIJ         *abt=(Mat_SeqAIJ*)submatBT->data;
      ai= abt->i;
      aj= abt->j;
      vals=abt->a;
      (*nzbs)[j]=abt->nz;
      (*bns)[j]=submatBT->rmap->n;
      (*bms)[j]=submatBT->cmap->n;
      //printf("j %d btn %d bms %d nz %d\n",j,submatBT->cmap->n,(*bms)[j],(*nzbs)[j]);
      for(j1=0; j1<(*bns)[j]; j1++) {
        (*irnbs)[j4]=ai[j1];
        j4++;
      }
      for(j1=0; j1<(*nzbs)[j]; j1++) {
        (*jcnbs)[j3]=aj[j1];
        //printf("jcn %d\n",(*jcnbs)[j3]);
        (*valbs)[j3]=vals[j1];
        j3++;
      }
      MatDestroy_SeqAIJ(submatBT);
    }
    free(submatB);
  }
  if(rank==0)for(j1=0; j1<j0; j1++)printf("nzzzzzz %d rank %d\n",(*nzbs)[j1],rank);
  return 0;
}*/

bool spline(ha_cgetype* y,ha_cgetype* x,ha_cgetype sx0,ha_cgetype sxn,int size,ha_cgetype* w,int laA){//sizeof(x)=size+1
  uvadd lasize;
  ha_cgetype *matval= (ha_cgetype *) calloc (6*4*size,sizeof(ha_cgetype));
  ha_cgetype *bval= (ha_cgetype *) calloc (4*size,sizeof(ha_cgetype));
  int *irn= (int *) calloc (6*4*size,sizeof(int));
  int *jcn= (int *) calloc (6*4*size,sizeof(int));
  if(size<3){
    printf("Too few points!!!\n");
    return 0;
  }
  int i,j,mindx=0,bindx=0,i1,i2;
  int *insize=(int *) calloc (4,sizeof(int));
  insize[0]=4*size;
  insize[1]=4*size;
  //printf("x0 %lf y0 %lf\n",x[0],y[0]);
  
  for(i=1;i<size+1;i++){//ai+bixi+cixi2+dixi3=yi//indx 1:n
    //printf("x0 %lf y0 %lf\n",x[i],y[i]);
    matval[mindx]=1;
    irn[mindx]=i;
    jcn[mindx]=i;
    mindx++;
    matval[mindx]=x[i];
    irn[mindx]=i;
    jcn[mindx]=i+size;
    mindx++;
    matval[mindx]=x[i]*x[i];
    irn[mindx]=i;
    jcn[mindx]=i+2*size;
    mindx++;
    matval[mindx]=x[i]*x[i]*x[i];
    irn[mindx]=i;
    jcn[mindx]=i+3*size;
    mindx++;
    bval[bindx]=y[i];
    bindx++;
  }
  for(i=size;i<2*size;i++){//aip+bipxi+cipxi2+dipxi3=yi//indx 0:n-1
    i1=i+1;
    i2=i-size;
    matval[mindx]=1;
    irn[mindx]=i1;
    jcn[mindx]=i1-size;
    mindx++;
    matval[mindx]=x[i2];
    irn[mindx]=i1;
    jcn[mindx]=i1;
    mindx++;
    matval[mindx]=x[i2]*x[i2];
    irn[mindx]=i1;
    jcn[mindx]=i1+size;
    mindx++;
    matval[mindx]=x[i2]*x[i2]*x[i2];
    irn[mindx]=i1;
    jcn[mindx]=i1+2*size;
    mindx++;
    bval[bindx]=y[i2];
    bindx++;
  }
    matval[mindx]=1;
    irn[mindx]=2*size+1;
    jcn[mindx]=1+size;
    mindx++;
    matval[mindx]=2*x[0];
    irn[mindx]=2*size+1;
    jcn[mindx]=1+2*size;
    mindx++;
    matval[mindx]=3*x[0]*x[0];
    irn[mindx]=2*size+1;
    jcn[mindx]=1+3*size;
    mindx++;
    bval[bindx]=(y[1]-y[0])/(x[1]-x[0]);//sx0;//
    bindx++;
  for(i=2*size+2;i<3*size+1;i++){//indx 1:n-1
    matval[mindx]=1;
    irn[mindx]=i;
    jcn[mindx]=i-size-1;
    mindx++;
    matval[mindx]=-1;
    irn[mindx]=i;
    jcn[mindx]=i-size;
    mindx++;

    matval[mindx]=2*x[i-2*size-1];
    irn[mindx]=i;
    jcn[mindx]=i-1;
    mindx++;
    matval[mindx]=-2*x[i-2*size-1];
    irn[mindx]=i;
    jcn[mindx]=i;
    mindx++;

    matval[mindx]=3*x[i-2*size-1]*x[i-2*size-1];
    irn[mindx]=i;
    jcn[mindx]=i+size-1;
    mindx++;
    matval[mindx]=-3*x[i-2*size-1]*x[i-2*size-1];
    irn[mindx]=i;
    jcn[mindx]=i+size;
    mindx++;
    bval[bindx]=0;
    bindx++;
  }
  for(i=3*size+1;i<4*size;i++){//indx 1:n-1
    matval[mindx]=2;
    irn[mindx]=i;
    jcn[mindx]=i-size;
    mindx++;
    matval[mindx]=-2;
    irn[mindx]=i;
    jcn[mindx]=i-size+1;
    mindx++;

    matval[mindx]=6*x[i-3*size];
    irn[mindx]=i;
    jcn[mindx]=i;
    mindx++;
    matval[mindx]=-6*x[i-3*size];
    irn[mindx]=i;
    jcn[mindx]=i+1;
    mindx++;
    bval[bindx]=0;
    bindx++;
  }
    matval[mindx]=1;
    irn[mindx]=4*size;
    jcn[mindx]=2*size;
    mindx++;
    matval[mindx]=2*x[size];
    irn[mindx]=4*size;
    jcn[mindx]=3*size;
    mindx++;
    matval[mindx]=3*x[size]*x[size];
    irn[mindx]=4*size;
    jcn[mindx]=4*size;
    mindx++;
    bval[bindx]=(y[size]-y[size-1])/(x[size]-x[size-1]);//sxn;//
    printf("slop %lf app slop %lf y1 %lf y0 %lf\n",sxn,bval[bindx],y[size],y[size-1]);
    bindx++;
  insize[2]=mindx;
  lasize=ceil((laA*10/100.0)*mindx);
  insize[3]=laA*10;
  irn=realloc(irn,lasize*sizeof(int));
  jcn=realloc(jcn,lasize*sizeof(int));
  matval=realloc(matval,lasize*sizeof(ha_cgetype));
  spec48_ssol2la_(insize,irn,jcn,matval,bval,w);
    printf("bindx %d val %lf w1 %lf w2 %lf w3 %lf x %lf\n",bindx,w[size-1]+w[size-1+size]*x[size]+w[size-1+2*size]*x[size]*x[size]+w[size-1+3*size]*x[size]*x[size]*x[size],w[size-1],w[size-1+size],w[size-1+2*size],y[size]);
  free(insize);
  free(matval);
  free(bval);
  free(irn);
  free(jcn);
  return 1;
}

bool Johansen(PetscBool nohsl,PetscInt VecSize,Mat A,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,Mat B,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec vecb,Vec vece,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar **ha_cofvar2,uvadd ncofvar,uvadd ncofele,uvadd nvarele,ha_cgeexovar **ha_cgeshock2,uvadd alltimeset,uvadd allregset,uvadd nintraeq,uvdim matsol,PetscInt Istart,PetscInt Iend,  uvadd nreg, uvadd ntime, uvadd *ha_eqadd, uvadd ndblock, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,uvdim laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol,uvdim nesteddbbd,int localsize,PetscInt *ndbbddrank1,forint* indata,uvdim mc66,forint *ptx,struct timeval begintime,ha_cgetype **xcf2){ //Johansen
  char tempfilenam[256],tempchar[256];
  PetscScalar value,*vals=NULL;
  PetscErrorCode ierr;
  PetscInt count,nz01,*ai=NULL,*aj=NULL;
  forint k=0,m=1;
  uvadd i,j,lasize;
  ha_cgetype *b1=NULL,*x0=NULL;
  bool IsIni;
  FILE* tempvar;
  ha_cgeexovar *ha_cgeshock;
  ha_cgeshock=*ha_cgeshock2;
  ha_cgetype *xcf;
  xcf=*xcf2;
  //printf("v0 %d p0 %d\n",ha_cgeshock[10].ExoIndx,(*ha_cgeshock2)[10].ExoIndx);

  //Vec vece,vecb;
  //PetscLogDouble time1,time0;
  clock_t timestr,timeend,timemulti;
  struct timeval endtime;
  struct timespec gettime_now,gettime_beg,gettime_end;
  long int start_time=0;
  double rep_time;
  size_t freadresult;
  ha_cgevar *ha_cofvar;
  ha_cofvar=*ha_cofvar2;
  ha_cgevar *ha_cofvar1=NULL;
  //printf("v0 %lf p0 %lf\n",ha_cofvar[10].varval,(*ha_cofvar2)[10].varval);
  
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
      free(*ha_cgeshock2);//
      //*ha_cgeshock2=(ha_cgeexovar*)realloc (*ha_cgeshock2,1*sizeof(ha_cgeexovar));
      *ha_cgeshock2=NULL;//realloc (ha_cgeshock,1*sizeof(ha_cgeexovar));
      ha_cgeshock=*ha_cgeshock2;
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
      free(*ha_cofvar2);//
      //*ha_cofvar2=(ha_cgevar*)realloc (*ha_cofvar2,1*sizeof(ha_cgevar));
      *ha_cofvar2=NULL;//realloc (ha_cofvar,1*sizeof(ha_cgevar));
      ha_cofvar=*ha_cofvar2;
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
    //ierr = PetscGetCPUTime(&time1);
    //CHKERRQ(ierr);
    //PetscViewerASCIIOpen(PETSC_COMM_WORLD, "B.txt", &viewer);
    //ierr = MatView(B,viewer);
    //CHKERRQ(ierr);
    //MatView(B,0);
    //printf("Vector e:\n");
    //VecView(e,0);
    ierr = VecDuplicate(vece,&vecb);
    CHKERRQ(ierr);
    ierr = MatMult(B,vece,vecb);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(vecb);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vecb);
    CHKERRQ(ierr);
    ierr = MatDestroy(&B);
    CHKERRQ(ierr);
    ierr = VecDestroy(&vece);
    CHKERRQ(ierr);
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
      clock_gettime(CLOCK_REALTIME, &gettime_beg);
      //ierr = PetscGetCPUTime(&time1);
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
      time(&timeend);
      //ierr = PetscGetCPUTime(&time1);
      //ierr = PetscPrintf(PETSC_COMM_WORLD,"One step solution %f\n",time1-time0);
      //if(rank==0)printf("One step calculation time %f\n",difftime(timeend,timestr));
      gettimeofday(&endtime, NULL);
      clock_gettime(CLOCK_REALTIME, &gettime_end);
      rep_time = ((double)(gettime_end.tv_nsec-gettime_beg.tv_nsec))/1000000000.0;
      if(rank==0)printf("One step calculation time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
      if(rank==0)printf("One step calculation real time %lf\n",rep_time);
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
      //ierr = PetscGetCPUTime(&time0);
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
        free(b1);//b1=realloc (b1,sizeof(ha_cgetype));
        b1=NULL;
      }
      else {
        x0=realloc (x0,VecSize*sizeof(ha_cgetype));
        //b1=realloc (b1,VecSize*sizeof(ha_cgetype));
        lasize=ceil((laA/100.0)*count);
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
        //b1=realloc (b1,sizeof(ha_cgetype));
        ierr = VecDestroy(&vecb);
        CHKERRQ(ierr);
      }
      //time(&timeend);
      //ierr = PetscGetCPUTime(&time1);
      //ierr = PetscPrintf(PETSC_COMM_WORLD,"One step solution %f\n",time1-time0);
      //if(rank==0)printf("One step calculation time %f\n",difftime(timeend,timestr));
      gettimeofday(&endtime, NULL);
      if(rank==0)printf("One step calculation time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
    }
    if(rank==rank_hsl) {
      if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
        printf("Error opening file\n");
      }
      *ha_cofvar2=(ha_cgevar*)realloc (*ha_cofvar2,(ncofele+nvarele)*sizeof(ha_cgevar));
      //ha_cofvar=(ha_cgevar*)calloc (ncofele+nvarele,sizeof(ha_cgevar));
      freadresult=fread(*ha_cofvar2, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
      fclose(tempvar);
      remove(tempfilenam);
      ha_cofvar=*ha_cofvar2;

      strcpy(tempfilenam,temdir);
      strcat(tempfilenam,"_tempshock");
      sprintf(tempchar, "%d",rank);
      strcat(tempfilenam,tempchar);
      strcat(tempfilenam,".bin");
      if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
        printf("Error opening file\n");
      }
      *ha_cgeshock2=(ha_cgeexovar*)realloc (*ha_cgeshock2,(nvarele)*sizeof(ha_cgeexovar));
      //ha_cgeshock=(ha_cgeexovar*)calloc (nvarele,sizeof(ha_cgeexovar));
      freadresult=fread(*ha_cgeshock2, sizeof(ha_cgeexovar),nvarele, tempvar);
      fclose(tempvar);
      remove(tempfilenam);
      ha_cgeshock=*ha_cgeshock2;      
    }
  //printf("v1 %d p1 %d\n",ha_cgeshock[10].ExoIndx,(*ha_cgeshock2)[10].ExoIndx);
  //printf("v1 %lf p1 %lf\n",ha_cofvar[10].varval,(*ha_cofvar2)[10].varval);
    //for(i=0; i<VecSize; i++)printf("x %f\n",x1[i]);
    //ha_cgetype *xs= (ha_cgetype *) calloc (nvarele+ncofele,sizeof(ha_cgetype));
    //ha_cgetype *xc= (ha_cgetype *) calloc (nvarele,sizeof(ha_cgetype));
    *xcf2=(ha_cgetype*)realloc (*xcf2,nvarele*sizeof(ha_cgetype));
    xcf=*xcf2;
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
      //printf("Rank %d Hello world1a!\n",rank);
      hnew_update(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele);
      strcpy(commsyntax,"formula");
      IsIni=false;
      //printf("Rank %d Hello world1b!\n",rank);
      hnew_calcff(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,IsIni);
      //printf("Rank %d Hello world1c!\n",rank);

      /*if(sol==0)for(i=0; i<ncofele+nvarele; i++) {
          xs[i]=ha_cofvar[i].varval;
        }*/
      //if(sol==0)
//      for(i=0; i<nvarele; i++) {
//        xcf[i]=varchange[i];
//      }
//      free(varchange);
    }
    ha_cofvar1=NULL;
    free(x0);
    //free(b1);
    //printf("Rank %d Hello world2!\n",rank);
    return 1;
  }

bool ModMidPoint(PetscBool nohsl,PetscInt VecSize,Mat* A1,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,Mat* B1,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec* vecb1,Vec *vece1,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar **ha_cofvar2,uvadd ncofvar,uvadd ncofele,uvadd nvarele,ha_cgeexovar **ha_cgeshock2,uvadd alltimeset,uvadd allregset,uvadd nintraeq,uvdim matsol,PetscInt Istart,PetscInt Iend,  uvadd nreg, uvadd ntime, uvadd *ha_eqadd, uvadd ndblock, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,uvdim laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol,uvdim nesteddbbd,int localsize,PetscInt *ndbbddrank1,forint* indata,uvdim mc66,forint *ptx,struct timeval begintime,uvdim subints,MPI_Fint fcomm,ha_cgetype **xcf2,int Isbiupd){ //Modified midpoint Pearson 1991
  char tempfilenam[256],tempchar[256],solchar[255];
  PetscScalar value,*vals;
  PetscErrorCode ierr;
  PetscInt count,nz01,*ai,*aj;
  forint k=0,m=1;
  forint tindx1;//,tindx2;
  ha_cgetype temp1,temp2;
  uvadd i,j,lasize;
  uvdim subindx;
  ha_cgetype *b1=NULL;
  ha_cgetype *x1=NULL;
  ha_cgetype *xcf;
  xcf=*xcf2;
  bool IsIni;
  FILE* tempvar;
  //Vec vece,vecb;
  PetscLogDouble time1,time0;
  clock_t timestr,timeend,timemulti;
  struct timeval endtime;
  struct timespec gettime_now,gettime_beg,gettime_end;
  long int start_time=0;
  double rep_time;
  size_t freadresult;
  ha_cgevar *ha_cofvar;
  ha_cofvar=*ha_cofvar2;
  ha_cgevar *ha_cofvar1;
  ha_cgeexovar *ha_cgeshock;
  ha_cgeshock=*ha_cgeshock2;
  Vec vece,vecb;
  Mat A,B;
  A=*A1;
  B=*B1;
  vece=*vece1,
  vecb=*vecb1;
  int stepcount;
  int nsteps=3;
  int sol;
  ha_cgetype vpercents=1.0,perprecis=0;
  FILE* solution;
  int maxsol=3;
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
    ha_cgetype *xc0= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
    //ha_cgetype *xc012= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
    //ha_cgetype *xc024= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
    ha_cgetype *xc12= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
    ha_cgetype *xc24= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
    int *xc124= (int *) calloc (1,sizeof(int));
    ha_cgetype *clag1= (ha_cgetype *) calloc (nvarele,sizeof(ha_cgetype));
    ha_cgetype *varchange= (ha_cgetype *) calloc (nvarele,sizeof(ha_cgetype));
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
          printf("rank %d subint %d sol %d stepcount %d nsteps %d\n",rank,subindx,sol,stepcount,nsteps);
          //printf("rank %d OKA!!!\n",rank);
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
            printf("rank %d OK!!!\n",rank);
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
          printf("OK!!!\n");
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
          //printf("OKA!!! dnnz %d\n",onnz[0]);

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
          printf("OKB!!!\n");

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
            free(*ha_cofvar2);
            //*ha_cofvar2=(ha_cgevar*)realloc (*ha_cofvar2,1*sizeof(ha_cgevar));
            *ha_cofvar2=NULL;
            ha_cofvar=*ha_cofvar2;

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
            free(*ha_cgeshock2);
            //*ha_cgeshock2=(ha_cgeexovar*)realloc (*ha_cgeshock2,1*sizeof(ha_cgeexovar));
            *ha_cgeshock2=NULL;
            ha_cgeshock=*ha_cgeshock2;
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
            time(&timestr);

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

            time(&timeend);
            //memcpy (ha_eqadd,ha_eqadd1,VecSize*sizeof(uvadd));
            //ha_eqadd1=realloc(ha_eqadd1,1*sizeof(uvadd));
            MPI_Barrier(PETSC_COMM_WORLD);
            ierr = PetscGetCPUTime(&time1);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"One step solution %f\n",time1-time0);
            if(rank==0)printf("One step calculation time %f\n",difftime(timeend,timestr));
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
              lasize=ceil((laA/100.0)*count);
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
            *ha_cgeshock2=(ha_cgeexovar*)realloc (*ha_cgeshock2,(nvarele)*sizeof(ha_cgeexovar));
            freadresult=fread(*ha_cgeshock2, sizeof(ha_cgeexovar),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            ha_cgeshock=*ha_cgeshock2;

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempvar");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            *ha_cofvar2=(ha_cgevar*)realloc (*ha_cofvar2,(ncofele+nvarele)*sizeof(ha_cgevar));
            freadresult=fread(*ha_cofvar2, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            ha_cofvar=*ha_cofvar2;

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
          printf("sol %d stepcount %d\n\n",sol,stepcount);
          MPI_Barrier(PETSC_COMM_WORLD);
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
          free(x1);
          x1=NULL;
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = VecAssemblyBegin(vece);
          CHKERRQ(ierr);
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = VecAssemblyEnd(vece);
          CHKERRQ(ierr);
          //x1=realloc (x1,sizeof(ha_cgetype));
          if(rank==rank_hsl) {
            if(stepcount==0) {
              hnew_update(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele);
            }
            else {
              hnew_mupdate(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele);
            }
            if(Isbiupd==1)hnew_biupd(rank,tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,nvarele,laA,subints,1,0,nsteps);
            if(Isbiupd==2)hnew_biupd(rank,tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,nvarele,laA,subints,1,2,nsteps);
            //printf("X2\n");
            strcpy(commsyntax,"formula");
            IsIni=false;
            hnew_calcff(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,IsIni);
            //printf("X3\n");
          }
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

        if(nohsl)MPI_Barrier(PETSC_COMM_WORLD);
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
          free(*ha_cofvar2);
          //*ha_cofvar2=(ha_cgevar*)realloc (*ha_cofvar2,1*sizeof(ha_cgevar));
          *ha_cofvar2=NULL;
          ha_cofvar=*ha_cofvar2;

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
          free(*ha_cgeshock2);
          //*ha_cgeshock2=(ha_cgeexovar*)realloc (*ha_cgeshock2,1*sizeof(ha_cgeexovar));
          *ha_cgeshock2=NULL;
          ha_cgeshock=*ha_cgeshock2;
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
          time(&timestr);

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

          time(&timeend);
          //memcpy (ha_eqadd,ha_eqadd1,VecSize*sizeof(uvadd));
          //ha_eqadd1=realloc(ha_eqadd1,1*sizeof(uvadd));
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time1);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"One step solution %f\n",time1-time0);
          if(rank==0)printf("One step calculation time %f\n",difftime(timeend,timestr));
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
            lasize=ceil((laA/100.0)*count);
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
          *ha_cgeshock2=(ha_cgeexovar*)realloc (*ha_cgeshock2,(nvarele)*sizeof(ha_cgeexovar));
          freadresult=fread(*ha_cgeshock2, sizeof(ha_cgeexovar),nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);
          ha_cgeshock=*ha_cgeshock2;

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempvar");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error opening file\n");
          }
          *ha_cofvar2=(ha_cgevar*)realloc (*ha_cofvar2,(ncofele+nvarele)*sizeof(ha_cgevar));
          freadresult=fread(*ha_cofvar2, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);
          ha_cofvar=*ha_cofvar2;

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
            *xcf2=(ha_cgetype*)realloc (*xcf2,(nvarele)*sizeof(ha_cgetype));
            xcf=*xcf2;
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
            *xcf2=(ha_cgetype*)realloc (*xcf2,nvarele*sizeof(ha_cgetype));
            xcf=*xcf2;
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
                    //xc24[k]=xcf[k];
                    //xcf[k]+=varchange[k]/45;
                    xc12[k]=xcf[k]-varchange[k]*kval1;
                    xc24[k]=xcf[k];
                    xcf[k]+=varchange[k]*kval2;
                    //}
                  }
                }
                else {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    //if(ha_cgeshock[ha_var[i].begadd+j].ShockId==0) {
                    //k=ha_var[i].begadd+j;
                    //if(k==1287)printf("sol!!!!!!!!!!!!!!!!!! %d step %d xc %lf xc0 %lf k %d\n",sol,stepcount,xc[k],xc0[k],k);
                    //xc24[k]=xcf[k];
                    //xc12[k]=xcf[k]-varchange[k]*xc0[k]/3;
                    //xcf[k]+=varchange[k]/45*xc0[k];
                    xc24[k]=xcf[k];
                    xc12[k]=xcf[k]-varchange[k]*xc0[k]*kval1;
                    xcf[k]+=varchange[k]*xc0[k]*kval2;//(100+xc0[k])*(100+varchange[k]/45)/100-100;//varchange[k]/45;
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
                    //xcf[k]-=20*varchange[k]/45*xc0[k];//
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
                    //xcf[k]+=64*varchange[k]/45*xc0[k];//
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
          if(Isbiupd==1)hnew_biupd(rank,tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,nvarele,laA,subints,1,0,nsteps);
          if(Isbiupd==2)hnew_biupd(rank,tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,nvarele,laA,subints,1,2,nsteps);
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
            free(*xcf2);
            *xcf2=NULL;
            xcf=*xcf2;
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
              //free(x0);
              ha_cofvar1=NULL;
              if(x1!=NULL)free(x1);
  }
