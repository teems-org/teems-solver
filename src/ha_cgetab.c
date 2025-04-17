//#include <ha_cgetab.h>
//#include <ha_cgefparse.h>
//#include "petscksp.h"
#include <ha_cgeglobal.h>

int ha_cgerecovar(char *fomulain) {
  int index,i,i1,j;
  char fpart1[TABREADLINE],*p=NULL;//,*p1=NULL;
  while (ha_cgefrstr(fomulain, " ", ""));
  p=strchr(fomulain,')');
  if (p==NULL) {
    while (ha_cgefrchr(fomulain, '[', '('));
    while (ha_cgefrchr(fomulain, ']', ')'));
    return 0;
  }
  strncpy(fpart1, fomulain, p-fomulain);
  fpart1[p-fomulain] = '\0';
  index=strrchr(fpart1,'(')-fpart1;//ha_cgerevfind(fpart1,"(");
  //if (index==0) {
  //fomulain[p-fomulain]=']';
  //fomulain[index]='[';
  //return 1;
  //}
  fpart1[index]='\0';
  for(i=index-1; i>-1; i--) {
    if(fpart1[i]=='+'||fpart1[i]=='-'||fpart1[i]=='*'||fpart1[i]=='/'||fpart1[i]=='['||fpart1[i]=='('||fpart1[i]==','||fpart1[i]=='^') {
      break;
    }
  }
  //p1=strpbrk(fpart1,"[*+/-(,^");
  if (i>-1) {
    i=index-i;
    switch (i) {
    case 1:
      fomulain[p-fomulain]=']';
      fomulain[index]='[';
      break;
    case 3:
      j=index-i;
      i=0;
      if(fpart1[j+1]=='i') {
        i++;
      }
      if(fpart1[j+2]=='f') {
        i++;
      }
      if (i==2) {
        fomulain[p-fomulain]=']';
        fomulain[index]='[';
      } else {
        fomulain[p-fomulain]='}';
        fomulain[index]='{';
      }
      break;
    case 4:
      j=index-i;
      i=0;
      if(fpart1[j+1]=='s') {
        i++;
      }
      if(fpart1[j+2]=='u') {
        i++;
      }
      if(fpart1[j+3]=='m') {
        i++;
      }
      i1=0;
      if(fpart1[j+1]=='a') {
        i1++;
      }
      if(fpart1[j+2]=='b') {
        i1++;
      }
      if(fpart1[j+3]=='s') {
        i1++;
      }
      if (i==3||i1==3) {
        fomulain[p-fomulain]=']';
        fomulain[index]='[';
      } else {
        fomulain[p-fomulain]='}';
        fomulain[index]='{';
      }
      break;
    case 5:
      j=index-i;
      i=0;
      if(fpart1[j+1]=='i') {
        i++;
      }
      if(fpart1[j+2]=='d') {
        i++;
      }
      if(fpart1[j+3]=='0') {
        i++;
      }
      if(fpart1[j+4]=='1') {
        i++;
      }
      i1=0;
      if(fpart1[j+1]=='l') {
        i1++;
      }
      if(fpart1[j+2]=='o') {
        i1++;
      }
      if(fpart1[j+3]=='g') {
        i1++;
      }
      if(fpart1[j+4]=='e') {
        i1++;
      }
      if (i==4||i1==4) {
        fomulain[p-fomulain]=']';
        fomulain[index]='[';
      } else {
        fomulain[p-fomulain]='}';
        fomulain[index]='{';
      }
      break;
    default:
      fomulain[p-fomulain]='}';
      fomulain[index]='{';
      break;
    }
  } else {
    i=index-i;
    switch (i) {
    case 1:
      fomulain[p-fomulain]=']';
      fomulain[index]='[';
      break;
    case 3:
      i=0;
      if(fpart1[0]=='i') {
        i++;
      }
      if(fpart1[1]=='f') {
        i++;
      }
      if (i==2) {
        fomulain[p-fomulain]=']';
        fomulain[index]='[';
      } else {
        fomulain[p-fomulain]='}';
        fomulain[index]='{';
      }
      break;
    case 4:
      i=0;
      if(fpart1[0]=='s') {
        i++;
      }
      if(fpart1[1]=='u') {
        i++;
      }
      if(fpart1[2]=='m') {
        i++;
      }
      i1=0;
      if(fpart1[0]=='a') {
        i1++;
      }
      if(fpart1[1]=='b') {
        i1++;
      }
      if(fpart1[2]=='s') {
        i1++;
      }
      if (i==3||i1==3) {
        fomulain[p-fomulain]=']';
        fomulain[index]='[';
      } else {
        fomulain[p-fomulain]='}';
        fomulain[index]='{';
      }
      break;
    case 5:
      i=0;
      if(fpart1[0]=='i') {
        i++;
      }
      if(fpart1[1]=='d') {
        i++;
      }
      if(fpart1[2]=='0') {
        i++;
      }
      if(fpart1[3]=='1') {
        i++;
      }
      if (i==4) {
        fomulain[p-fomulain]=']';
        fomulain[index]='[';
      } else {
        fomulain[p-fomulain]='}';
        fomulain[index]='{';
      }
      break;
    default:
      fomulain[p-fomulain]='}';
      fomulain[index]='{';
      break;
    }
  }
  return 1;
}

int ha_cgecutsum(char *formula) {
  char *t1=")",*t2="(";
  int i,l,count=0;
  l=strchr(formula,'(')-formula;//ha_cgefind(formula,"(");

  for (i=l; formula[i]; i++) {
    if (formula[i]==*t1) {
      count=count-1;
    }
    if (formula[i]==*t2) {
      count=count+1;
    }
    //printf("i %d count %d",i,count);
    if (count==0) {
      break;
    }
  }
  formula[i+1]='\0';
  return 1;
}
/*
int hcge_matrele(char *fname, char *commsyntax,ha_cgeset *ha_set,int nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof, long int ncof,hcge_cof *ha_var, long int nvar, ha_cgevar *ha_cofvar, long int ncofvar, long int ncofele,ha_cgeexovar *ha_cgeshock,PetscInt Istart,PetscInt Iend,PetscInt *dnz,PetscInt *dnnz,PetscInt *onz,PetscInt *onnz,PetscInt *dnzB,PetscInt *dnnzB,PetscInt *onzB,PetscInt *onnzB)
{
  FILE * filehandle;
  char line[TABREADLINE],line1[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[NAMESIZE];//,var1[NAMESIZE];
  unsigned short int i,i1=0,i2=0,l,l1=0,fdim;//,totalsum,sumcount=1,c1,c2,svar1,svar2;
  PetscInt index;//j=0,varadd,
  char *readitem=NULL,*p=NULL,*p1=NULL,*p2=NULL;//,*p3=NULL,*p4=NULL;
  long int dcount,dcountdim1[MAXVARDIM],nloops,i3;//,nsumele,l2;
  unsigned int np,nlinvars,sizelinvars,templinvars;//i4,
  //int sumindx;
  //char varname[NAMESIZE];
  uvdim i4;

  filehandle = fopen(fname,"r");
  index=0;

  while (ha_cgertabl(commsyntax,filehandle,line)) {
    if (strstr(line,"(DEFAULT")==NULL) {
      ha_cgefrstr1(line, commsyntax, "");
      while (ha_cgefrstr(line," ", ""));
      while (ha_cgefrchr(line, '[', '('));
      while (ha_cgefrchr(line, ']', ')'));
      while (ha_cgefrchr(line, '{', '('));
      while (ha_cgefrchr(line, '}', ')'));
      strcpy(linecopy,line);
      //printf("equation %s\n",line);
      fdim=ha_cgenfind(line, "(all,");
      ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (fdim,sizeof(ha_cgesetindx));
      i=ha_cgerevfind(line, "(all,");
      readitem=&line[i];
      readitem = strtok(readitem,",");
      readitem = strtok(NULL,")");
      readitem = strtok(NULL,"=");
      strcpy(vname,readitem);
      strcpy(line,linecopy);
      readitem = strtok(line,"=");
      readitem = strtok(NULL,";");
      strcat(readitem,"-");
      strcat(readitem,vname);
      while (ha_cgerecovar(readitem)==1);
      strcpy(line1,readitem);
      np=ha_cgenfind(readitem,"p_");
      hcge_cof *LinVars= (hcge_cof *) calloc (np,sizeof(hcge_cof));
      unsigned long int *LinVarsBA= (unsigned long int *) calloc (np,sizeof(unsigned long int));
      i3=0;
      sizelinvars=0;
      for (i=0; i<np; i++) {
        readitem=strstr(readitem,"p_");
        readitem++;
        readitem++;
        p=strpbrk(readitem,"{+*-/^)");
        if(p!=NULL) {
          strncpy(vname,readitem,p-readitem);
          vname[p-readitem]='\0';
          i2=0;
          for (i1=0; i1<i+1; i1++) if (strcmp(LinVars[i1].cofname,vname)==0) {
              i2=1;
              break;
            }
          if (i2==0) {
            strcpy(LinVars[i3].cofname,vname);
            for (l=0; l<nvar; l++) {
              if (strcmp(ha_var[l].cofname,vname)==0) {
                LinVars[i3].begadd=ha_var[l].begadd;
                LinVars[i3].size=ha_var[l].size;
                templinvars=1;
                for(l1=0; l1<ha_var[l].size; l1++) {
                  strcpy(LinVars[i3].dimsets[l1],ha_var[i3].dimsets[l1]);
                  LinVars[i3].dims[l1]=ha_var[l].dims[l1];
                  templinvars=templinvars*ha_var[l].dims[l1];
                }
                LinVarsBA[i]=templinvars;
                sizelinvars=sizelinvars+templinvars;
                break;
              }
            }
            i3++;
          }
          //LinVars[i][p-readitem]='\0';
        } else {
          strcpy(vname,readitem);
          for (i1=0; i1<i+1; i1++) if (strcmp(LinVars[i1].cofname,vname)==0) {
              i2=1;
              break;
            }
          if (i2==0) {
            strcpy(LinVars[i3].cofname,vname);
            for (l=0; l<nvar; l++) {
              if (strcmp(ha_var[l].cofname,vname)==0) {
                LinVars[i3].begadd=ha_var[l].begadd;
                LinVars[i3].size=ha_var[l].size;
                templinvars=1;
                for(l1=0; l1<ha_var[l].size; l1++) {
                  strcpy(LinVars[i3].dimsets[l1],ha_var[i3].dimsets[l1]);
                  LinVars[i3].dims[l1]=ha_var[l].dims[l1];
                  templinvars=templinvars*ha_var[l].dims[l1];
                }
                LinVarsBA[i]=templinvars;
                sizelinvars=sizelinvars+templinvars;
                break;
              }
            }
            i3++;
          }
          break;
        }
      }
      nlinvars=i3;
      //for (i=0; i<np; i++) {
      //printf("lin %s size %d\n",LinVars[i].cofname,LinVarsBA[i]);
      //}

      //printf("readitem %s\n",readitem);
      strcpy(line,linecopy);

      nloops=1;
      if (fdim>0) {
        for (i=0; i<fdim; i++) {
          if(i==0) {
            readitem = strtok(line,",");
          } else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
          strcpy(arSet[i].arIndx,readitem);
          readitem = strtok(NULL,")");
          strcpy(arSet[i].arSet,readitem);
          for (i4=0; i4<nset; i4++) if(strcmp(readitem,ha_set[i4].setname)==0) {
              arSet[i].SetSize=ha_set[i4].size;
              break;
            }
          nloops=nloops*ha_set[i4].size;
          dcount=fdim-i;
          if(dcount==fdim) {
            dcountdim1[dcount-1]=1;
          } else {
            dcountdim1[dcount-1]=ha_set[i4].size*dcountdim1[dcount];
          }
        }
      }
      *dnz=0;
      *onz=0;
      //printf("lin1 %s\n",line1);
      i1=0;
      for (l1=0; l1<nlinvars; l1++) {
        strcpy(line,&line1[i1]);
        //printf("lin1 %s\n",line);
        readitem=strstr(line,LinVars[l1].cofname);
        //printf("lin1 %s\n",readitem);
        if (LinVars[l1].size>0) {
          for (i=0; i<LinVars[l1].size-1; i++) {
            strcpy(p,"sum(");
            readitem=strchr(readitem,'{');
            readitem++;
            p2=strchr(readitem,',');
            p2++;
            //printf("lin1 %s\n",p2);
            strcpy(p1,readitem);
            //printf("p1 %s\n",p1);
            p1[p2-readitem]='\0';
            strcat(p,p1);
            //printf("lin1 %s\n",p);
            //p=strstr(line,p);
            //if (p!=NULL) printf("in sum %d\n",hcge_isinsum(line,readitem-line,p-line));
          }
          strcpy(p,"sum(");
          p1=strchr(readitem,'{');
          p1++;
          p2=strchr(p1,'}');
          p2++;
          p1[p2-p1]='\0';
          strcat(p,p1);
          //printf("lin1 %s\n",p);
          //p=strstr(line,p);
          //if (p!=NULL) printf("r-l %d p-l %d",readitem-line,p-line);
          //if (p!=NULL) printf("in sum %d\n",hcge_isinsum(line,readitem-line,p-line));
        }
        //ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (fdim,sizeof(ha_cgesetindx));


        //ha_cgesetindx *arSum= (ha_cgesetindx *) calloc (fdim,sizeof(ha_cgesetindx));



        for (l=0; l<nloops; l++) {
          if(Istart<=index&&index<Iend) {
          }
        }
        index=index+1;
      }
      free(LinVars);
      free(LinVarsBA);
      free(arSet);
    }
  }
  return 1;
}*/
/*
int hcge_isinsum(char *formulain, int addLin, int addSum)
{
  char *t1=")",*t2="(",*p;
  int i,l,count=0;
  p=&formulain[addSum];
  l=strchr(p,'(')-p;//ha_cgefind(formula,"(");

  for (i=l; p[i]; i++) {
    if (p[i]==*t1) {
      count=count-1;
    }
    if (p[i]==*t2) {
      count=count+1;
    }
    //printf("i %d count %d",i,count);
    if (count==0) {
      break;
    }
  }
  if (addLin>addSum+i+l) {
    return 0;
  } else {
    return 1;
  }
}*/

int hcge_repllin(char *formulain,int linindx) {
  char *p,*p1,*p2,leftline[TABREADLINE];//,*p3
  int np,i,l,pcheck;//,pl,mn;//d,
  np=ha_cgenfind(formulain,"p_");
  p=&formulain[0];
  p2=&formulain[0];
  for (i=0; i<np; i++) {
    //l=strlen(p);
    p=p+ha_cgefind(p,"p_");
    if(p==&formulain[0]) {
      pcheck=0;
      p1=p;
      while (pcheck==0) {
        p1=strpbrk(p1,"+-");
        if (p1==NULL) {
          strcpy(leftline,p2);
        } else {
          strncpy(leftline,p2,p1-p2);
          leftline[p1-p2]='\0';
        }
        if(ha_cgenfind(leftline,"(")==ha_cgenfind(leftline,")")) {
          pcheck=1;
        }
        p1++;
      }
      pcheck=0;
      p1=&leftline[0];
      p2=p1;
      l=0;
      while (p1!=NULL) {
        p1=strpbrk(p1,"+-");
        if (p1==NULL) {
          break;
        }
        if (ha_cgenfind(leftline,"(")==0&&ha_cgenfind(leftline,")")==0) {
          p2=p1;
          p1++;
          l++;
        } else {
          if (ha_cgenfind(p1,"(")!=ha_cgenfind(p1,")")) {
            l++;
            p1++;
          } else {
            p2=p1;
            p1++;
            l++;
          }
        }
      }
      break;
    } else {
      if(*(p-1)=='*'||*(p-1)=='+'||*(p-1)=='-'||*(p-1)=='('||*(p-1)==',') {
        pcheck=0;
        p1=p;
        while (pcheck==0) {
          p1=strpbrk(p1,"+-");
          if (p1==NULL) {
            strcpy(leftline,p2);
          } else {
            strncpy(leftline,p2,p1-p2);
            leftline[p1-p2]='\0';
          }
          if(ha_cgenfind(leftline,"(")==ha_cgenfind(leftline,")")) {
            pcheck=1;
          }
          p1++;
        }
        pcheck=0;
        p1=&leftline[0];
        p2=p1;
        l=0;
        while (p1!=NULL) {
          p1=strpbrk(p1,"+-");
          if (p1==NULL) {
            break;
          }
          if (ha_cgenfind(leftline,"(")==0&&ha_cgenfind(leftline,")")==0) {
            p2=p1;
            p1++;
            l++;
          } else {
            if (ha_cgenfind(p1,"(")!=ha_cgenfind(p1,")")) {
              l++;
              p1++;
            } else {
              p2=p1;
              p1++;
              l++;
            }
          }
        }
        break;
      }
    }
    p++;
  }
  //printf("p2 %s\n",p2);
  //if (l>0) p2--;
  //if (*p2=='+'||*p2=='-') {strcpy(formulain,"0");strcat(formulain,p2);}
  //else strcpy(formulain,p2);
  strcpy(formulain,p2);
  return 1;
}

int hcge_rlinzero(char *formulain,int linindx) {
  char *p,*p1=NULL;//,line1[TABREADLINE];
  int np,i,l=linindx;//,l,d;
  np=ha_cgenfind(formulain,"p_");
  p=&formulain[0];
  for (i=0; i<np; i++) {
    //l=strlen(p);
    p=p+ha_cgefind(p,"p_");
    //l=p-formulain;
    //p2=p--;
    //printf("p1 %s\n",p);

    if(p==&formulain[0]) {
      p1=strpbrk(p,"}+*-/^)");
      if (p1==NULL) {
        p1=strchr(p,'\0');
      } else if (*p1=='}') {
        p1++;
      }
      if (i==l) {
        p=p1;
        //printf("lin %d p %s\n",linindx,p);
      } else {
        //printf("for %s\n",formulain);
        if (*p1=='\0') {
          *p='0';
          p++;
          *p='\0';
        } else {
          memmove(p+1,p1,strlen(p1)+1);  //formulain[l-d+1]='\0';
          *p='0';
        }
      }
    } else {
      if(*(p-1)=='*'||*(p-1)=='+'||*(p-1)=='-'||*(p-1)=='('||*(p-1)==',') {
        p1=strpbrk(p,"}+*-/^)");
        if (p1==NULL) {
          p1=strchr(p,'\0');
        } else if (*p1=='}') {
          p1++;
        }
        if (i==l) {
          p=p1;
          //printf("lin %d p %s\n",linindx,p);
        } else {
          //printf("for %s\n",formulain);
          if (*p1=='\0') {
            *p++='0';
            p++;
            *p='\0';
          } else {
            memmove(p+1,p1,strlen(p1)+1);  //formulain[l-d+1]='\0';
            *p++='0';
          }
        }
      } else {
        p++;
        l++;
      }
    }
  }
  return 1;
}

uvadd hcge_readff(char *fname, int niodata, hcge_iodata *iodata, char *commsyntax,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele,hcge_cof *ha_cof,uvadd ncof, ha_cgecofele *ha_cofele,uvadd ncofele,hcge_cof *ha_var,uvadd nvar, ha_cgecofele *ha_varele,uvadd nvarele) {
  FILE * filehandle, * filehandle1;
  char line[DATREADLINE]="\0",linecopy[DATREADLINE],line1[TABREADLINE],*p,*p1;//,line1[DATREADLINE]
  char vname[NAMESIZE],setindx[NAMESIZE],header[NAMESIZE],*vname1=NULL,argu[TABREADLINE],varset[MAXVARDIM][NAMESIZE];//,varindx[MAXVARDIM][NAMESIZE],varindx1[MAXVARDIM][NAMESIZE];
  //char *commsyntax="SET";
  //unsigned short int l1,l2,l3,l4;//,vsize
  uvadd j=0,i,l1,count1,recount,recount1,n,n1,n2,dims,antidim[MAXVARDIM],dim[MAXVARDIM],begadd[MAXVARDIM],supsetid[MAXVARDIM],index[MAXVARDIM];
  int count2;//0 no finding;1 find var; 2 find header cof; 4 find header var;
  //long int antidim1[MAXVARDIM],dim1[MAXVARDIM],begadd1[MAXVARDIM],index1[MAXVARDIM],varsize[MAXVARDIM];
  char *readitem=NULL,*copyvar;
  ha_cgetype val;
  int k0,k1;

  filehandle = fopen(fname,"r");
  if(filehandle==NULL){
    printf("Error!!! No such %s file!\n",fname);
    return -1;
  }
  //fileout= fopen("my.txt", "wb");

  while (ha_cgertabl(commsyntax,filehandle,line,DATREADLINE)) {
    //printf("%s\n",line);
    /*while (ha_cgefrstr(line,"\n", ""));
    while (ha_cgefrstr(line,"\r", ""));
    while (ha_cgefrstr(line,"  ", " "));
    strcat(line, "\n");
    while (ha_cgedrcmt(line,"!"));
    while (ha_cgedrcmt(line,"#"));
    while (ha_cgefrstr(line,"  ", " "));*/
    strcpy(linecopy,line);
    //printf("line %s\n",line);
    k0=ha_cgefind(line,"from file ");
    if(k0>-1) {
      k1=ha_cgefind(line+k0+10," ");
      strncpy(line1,line+k0+10,k1);
      line1[k1]='\0';
      //printf("line1 %s\n",line1);
      for (k0=0; k0<niodata; k0++) if (strcmp(line1,iodata[k0].logname)==0) {
          break;
        }
      if(k0==niodata){
        printf("Error!!!! No such file %s in cmf file!\n",line1);
        return -1;
      }
    }
    if (strchr(line,'(')==NULL) {
      readitem = strtok(line," ");
      readitem = strtok(NULL," ");
      strcpy(vname,readitem);
      readitem = strtok(NULL,"\"");
      readitem = strtok(NULL,"\"");
      strcpy(header,readitem);
      n1=0;
      while (header[n1]!='\0'){
        if(header[n1]==' '){
          header[n1]='\0';
          break;
        }
        n1++;
      }
      //vsize=0;
      //while (vname[vsize] != '\0') vsize++;
      count2=0;
      for (i=0; i<nvar; i++) {
        vname1= strtok(ha_var[i].cofname,"(");
        //printf("vname %s varname %s varsize %d\n",vname,vname1,ha_var[i].size);
        if (strcmp(vname1,vname)==0) { //,vsize ha_var[i].cofname
          //printf("vname %s varname %s varsize %d\n",vname,vname1,ha_var[i].size);
          ha_var[i].suplval=true;
          dims=1;
          for (n1=0; n1<ha_var[i].size; n1++) {
            dims=dims*ha_set[ha_var[i].setid[n1]].size;
          }
          antidim[0]=ha_set[ha_var[i].setid[1]].size;
          antidim[1]=1;
          antidim[2]=antidim[0]*ha_set[ha_var[i].setid[0]].size;
          if (ha_var[i].size>3) {
            for (n1=3; n1<ha_var[i].size; n1++) {
              antidim[n1]=antidim[n1-1]*ha_set[ha_var[i].setid[n1-1]].size;
            }
          }
          strcpy(vname,readitem);
          n1=0;
          while (vname[n1]!='\0'){
            if(vname[n1]==' '){
              vname[n1]='\0';
              break;
            }
            n1++;
          }
          filehandle1 = fopen(iodata[k0].filname,"r");
          if(filehandle1==NULL){
            printf("Error!!! No such %s file!\n",iodata[k0].filname);
            return -1;
          }
          while (fgets(line,DATREADLINE,filehandle1)) {
            readitem = strtok(line,"\"");
            readitem = strtok(NULL,"\"");
            //printf("read %s\n",readitem);
            if (readitem != NULL) {
              //printf("read %s\n",readitem);
              //printf("vname %s\n",vname);
              if (strcmp(header,vname) == 0) {
                //printf("dim %s\n",readitem);
                count2=4;
                recount=0;
                while (fgets(line,DATREADLINE,filehandle1)) {
                  count1=0;
                  copyvar=line;
                  while (*copyvar==' ' || *copyvar=='\t') {
                    copyvar++;  // skip whitespaces
                  }
                  if (*copyvar=='\r') {
                    copyvar++;
                  }
                  if (*copyvar=='\n') {
                    count1=1;
                  }
                  if (count1!=1) {
                    n=ha_cgenchf(line,',');
                    if (n==0) {
                      readitem = strtok(line,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if(ha_var[i].size>1) {
                        for (n1=ha_var[i].size-1; n1>1; n1--) {
                          index[n1]=(uvadd) l1/antidim[n1];
                          //printf("indx %d\n",index[n1]);
                          l1=l1-antidim[n1]*index[n1];
                        }
                        index[0]=(uvadd) l1/antidim[0];
                        l1=l1-antidim[0]*index[0];
                        //printf("reacoun1 %d\n",recount1);
                        index[1]=(uvadd) l1/antidim[1];
                        l1=l1-antidim[1]*index[1];
                      } else {
                        index[0]=recount;
                      }
                      recount1=0;
                      for (n1=0; n1<ha_var[i].size; n1++) {
                        recount1=recount1+index[n1]*ha_var[i].antidims[n1];
                        //printf("reacoun %d ssh %d\n",recount1,ha_setele[begadd[n1]+index[n1]].setsh);
                      }
                      //printf("reacoun1 %d\n",recount1);
                      ha_varele[ha_var[i].begadd+recount1].cofval=val;
                      recount++;
                      if (recount>=dims) {
                        break;
                      }
                    } else {
                      for (n2=0; n2<n; n2++) {
                        //printf("n1 %d\n",n2);
                        if (n2==0) {
                          readitem = strtok(line,",");
                        } else {
                          readitem = strtok(NULL,",");
                        }
                        val=atof(readitem);//sscanf(readitem, "%lf", &val);
                        l1=recount;
                        if(ha_var[i].size>1) {
                          for (n1=ha_var[i].size-1; n1>1; n1--) {
                            index[n1]=(uvadd) l1/antidim[n1];
                            //printf("indx %d\n",index[n1]);
                            l1=l1-antidim[n1]*index[n1];
                          }
                          index[0]=(uvadd) l1/antidim[0];
                          l1=l1-antidim[0]*index[0];
                          index[1]=(uvadd) l1/antidim[1];
                          l1=l1-antidim[1]*index[1];
                          //printf("indx %d\n",index[1]);
                          //printf("indx %d\n",index[0]);
                        } else {
                          index[0]=recount;
                        }
                        recount1=0;
                        for (n1=0; n1<ha_var[i].size; n1++) {
                          recount1=recount1+index[n1]*ha_var[i].antidims[n1];
                          //printf("reacoun %d ssh %d\n",recount1,ha_setele[begadd[n1]+index[n1]].setsh);
                        }
                        //printf("reacoun1 %d\n",recount1);
                        ha_varele[ha_var[i].begadd+recount1].cofval=val;
                        recount++;
                      }
                      readitem = strtok(NULL,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if(ha_var[i].size>1) {
                        for (n1=ha_var[i].size-1; n1>1; n1--) {
                          index[n1]=(uvadd) l1/antidim[n1];
                          //printf("indx %d\n",index[n1]);
                          l1=l1-antidim[n1]*index[n1];
                        }
                        index[0]=(uvadd) l1/antidim[0];
                        l1=l1-antidim[0]*index[0];
                        index[1]=(uvadd) l1/antidim[1];
                        l1=l1-antidim[1]*index[1];
                      } else {
                        index[0]=recount;
                      }
                      recount1=0;
                      for (n1=0; n1<ha_var[i].size; n1++) {
                        recount1=recount1+index[n1]*ha_var[i].antidims[n1];
                        //printf("reacoun %d ssh %d\n",recount1,ha_setele[begadd[n1]+index[n1]].setsh);
                      }
                      //printf("reacoun1 %d\n",recount1);
                      ha_varele[ha_var[i].begadd+recount1].cofval=val;
                      //newrecount=rdatconv(ha_var,i,recount);
                      //printf("newrec %d rec %d\n",newrecount,recount);
                      recount++;
                      if (recount>=dims) {
                        break;
                      }
                    }
                  }

                }
                break;
              }

            }
            if(feof(filehandle1)) {
              break;
            }
          }
          if(count2!=4){
            printf("Error!!!! %s is not a header in data file!\n",header);
            return -1;
          }
          fclose(filehandle1);
          break;
        }
      }
      if(i<nvar)count2=1;
      for (i=0; i<ncof; i++) {
        vname1= strtok(ha_cof[i].cofname,"(");
        //printf("vname %s varname %s varsize %d\n",vname,vname1,ha_cof[i].size);
        if (strcmp(vname1,vname)==0) { //,vsize ha_cof[i].cofname
          //printf("i %ld vname1 %s varname %s varsize %d\n",i,vname,vname1,ha_cof[i].size);
          ha_cof[i].suplval=true;
          dims=1;
          for (n1=0; n1<ha_cof[i].size; n1++) {
            dims=dims*ha_set[ha_cof[i].setid[n1]].size;
          }
          antidim[0]=ha_set[ha_cof[i].setid[1]].size;
          antidim[1]=1;
          antidim[2]=antidim[0]*ha_set[ha_cof[i].setid[0]].size;
          if (ha_cof[i].size>3) {
            for (n1=3; n1<ha_cof[i].size; n1++) {
              antidim[n1]=antidim[n1-1]*ha_set[ha_cof[i].setid[n1-1]].size;
            }
          }
          //strcpy(vname,readitem);
//           n1=0;
//           while (vname[n1]!='\0'){
//             if(vname[n1]==' '){
//               vname[n1]='\0';
//               break;
//             }
//           }
          filehandle1 = fopen(iodata[k0].filname,"r");
          if(filehandle1==NULL){
            printf("Error!!! No such %s file!\n",iodata[k0].filname);
            return -1;
          }
          while (fgets(line,DATREADLINE,filehandle1)) {
            readitem = strtok(line,"\"");
            readitem = strtok(NULL,"\"");
            if (readitem != NULL) {
              strcpy(vname,readitem);
              //printf("read %s vname %s head %s strlen %ld\n",readitem,vname,header,strlen(vname));
              n1=0;
              while (vname[n1]!='\0'){
                if(vname[n1]==' '){
                  vname[n1]='\0';
                  break;
                }
                n1++;
              }
              //printf("read %s\n",readitem);
              //printf("vname %s\n",vname);
              if (strncmp(vname,header,strlen(header))==0&&strlen(header)==strlen(vname)) {
                //if(i==152)printf("dim %s vname %s head %s begadd %ld\n",readitem,vname,header,ha_cof[i].begadd);
                count2=2;
                recount=0;
                while (fgets(line,DATREADLINE,filehandle1)) {
                  //if(i==152)printf("line %s\n",line);
                  count1=0;
                  copyvar=line;
                  while (*copyvar==' ' || *copyvar=='\t') {
                    copyvar++;  // skip whitespaces
                  }
                  if (*copyvar=='\r') {
                    copyvar++;
                  }
                  if (*copyvar=='\n') {
                    count1=1;
                  }
                  if (count1!=1) {
                    n=ha_cgenchf(line,',');
                    if (n==0) {
                      readitem = strtok(line,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if(ha_cof[i].size>1) {
                        for (n1=ha_cof[i].size-1; n1>1; n1--) {
                          index[n1]=(uvadd) l1/antidim[n1];
                          //printf("indx %d\n",index[n1]);
                          l1=l1-antidim[n1]*index[n1];
                        }
                        index[0]=(uvadd) l1/antidim[0];
                        l1=l1-antidim[0]*index[0];
                        //printf("reacoun1 %d\n",recount1);
                        index[1]=(uvadd) l1/antidim[1];
                        l1=l1-antidim[1]*index[1];
                      } else {
                        index[0]=recount;
                      }
                      recount1=0;
                      for (n1=0; n1<ha_cof[i].size; n1++) {
                        recount1=recount1+index[n1]*ha_cof[i].antidims[n1];
                        //printf("reacoun %d ssh %d\n",recount1,ha_setele[begadd[n1]+index[n1]].setsh);
                      }
                      //printf("reacoun1 %d\n",recount1);
                      ha_cofele[ha_cof[i].begadd+recount1].cofval=val;
                      recount++;
                      if (recount>=dims) {
                        break;
                      }
                    } else {
                      for (n2=0; n2<n; n2++) {
                        if (n2==0) {
                          readitem = strtok(line,",");
                        } else {
                          readitem = strtok(NULL,",");
                        }
                        val=atof(readitem);//sscanf(readitem, "%lf", &val);
                        //printf("vname %s val %f\n",vname,val);
                        l1=recount;
                        if(ha_cof[i].size>1) {
                          for (n1=ha_cof[i].size-1; n1>1; n1--) {
                            index[n1]=(uvadd) l1/antidim[n1];
                            //printf("indx %d\n",index[n1]);
                            l1=l1-antidim[n1]*index[n1];
                          }
                          index[0]=(uvadd) l1/antidim[0];
                          l1=l1-antidim[0]*index[0];
                          index[1]=(uvadd) l1/antidim[1];
                          l1=l1-antidim[1]*index[1];
                          //printf("indx %d\n",index[1]);
                          //printf("indx %d\n",index[0]);
                        } else {
                          index[0]=recount;
                        }
                        recount1=0;
                        for (n1=0; n1<ha_cof[i].size; n1++) {
                          recount1=recount1+index[n1]*ha_cof[i].antidims[n1];
                          //printf("reacoun %d ssh %d\n",recount1,ha_setele[begadd[n1]+index[n1]].setsh);
                        }
                        //printf("reacoun1 %d\n",recount1);
                        ha_cofele[ha_cof[i].begadd+recount1].cofval=val;
                        recount++;
                      }
                      readitem = strtok(NULL,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if(ha_cof[i].size>1) {
                        for (n1=ha_cof[i].size-1; n1>1; n1--) {
                          index[n1]=(uvadd) l1/antidim[n1];
                          //printf("indx %d\n",index[n1]);
                          l1=l1-antidim[n1]*index[n1];
                        }
                        index[0]=(uvadd) l1/antidim[0];
                        l1=l1-antidim[0]*index[0];
                        index[1]=(uvadd) l1/antidim[1];
                        l1=l1-antidim[1]*index[1];
                      } else {
                        index[0]=recount;
                      }
                      recount1=0;
                      for (n1=0; n1<ha_cof[i].size; n1++) {
                        recount1=recount1+index[n1]*ha_cof[i].antidims[n1];
                        //printf("reacoun %d ssh %d\n",recount1,ha_setele[begadd[n1]+index[n1]].setsh);
                      }
                      //printf("reacoun1 %d\n",recount1);
                      ha_cofele[ha_cof[i].begadd+recount1].cofval=val;
                      //newrecount=rdatconv(ha_cof,i,recount);
                      //printf("newrec %d rec %d\n",newrecount,recount);
                      recount++;
                      if (recount>=dims) {
                        break;
                      }
                    }
                  }

                }
                break;
              }

            }
            if(feof(filehandle1)) {
              break;
            }
          }
          if(count2!=2){
            printf("Error1!!!! %s is not a header in data file!\n",header);
            return -1;
          }
          fclose(filehandle1);
          break;
        }
      }
      if(count2==0&&i==ncof){
        printf("Error!!!! %s is not a variable coefficient or parameter!\n",vname);
        return -1;
      }
    } else {
      n=ha_cgenchf(line,'(');
      readitem = strtok(line,")");
      for (n1=2; n1<n; n1++) {
        readitem = strtok(NULL,")");
      }
      //printf("line %s n %d\n",readitem,n);
      //readitem = strtok(NULL," ");
      readitem = strtok(NULL,"(");
      //printf("line %s\n",readitem);
      strcpy(vname,readitem);
      ha_cgefrstr(vname," ","");
      //printf("var %s argu %s\n",vname,argu);
      readitem = strtok(NULL,")");
      strcpy(argu,readitem);
      strcat(argu,",");
      readitem = strtok(NULL,"\"");
      readitem = strtok(NULL,"\"");
      strcpy(header,readitem);
      n1=0;
      while (header[n1]!='\0'){
        if(header[n1]==' '){
          header[n1]='\0';
          break;
        }
        n1++;
      }
      //printf("var %s argu %s header %s\n",vname,argu,header);
      strcpy(line,linecopy);
      n=ha_cgenchf(argu,',');
      for (n1=0; n1<n; n1++) {
        if (n1==0) {
          readitem = strtok(argu,",");
        } else {
          readitem = strtok(NULL,",");
        }
        strcpy(setindx,",");
        strcat(setindx,readitem);
        strcat(setindx,",");
        //printf("line %s\n",setindx);
        p=strstr(line,setindx);
        p1=strchr(p+1,',');
        p=strchr(p1,')');
        strncpy(varset[n1],p1+1,p-p1-1);
        varset[n1][p-p1-1]='\0';
        //strcpy(line,linecopy);
        //printf("setname %s\n",varset[n1]);
      }
      //printf("here!!!\n");
      //vsize=0;
      //while (vname[vsize] != '\0') vsize++;
      count2=0;
      for (i=0; i<ncof; i++) {
        vname1= strtok(ha_cof[i].cofname,"(");
        //printf("vname %s varname %s varsize %d\n",vname,vname1,ha_var[i].size);
        if (strcmp(vname1,vname)==0) { //,vsize ha_var[i].cofname
          dims=1;
          for (n1=0; n1<ha_cof[i].size; n1++) {
            /*for (n=0; n<nset; n++) {
              if (strcmp(varset[n1],ha_set[n].setname)==0) {
                dims=dims*ha_set[n].size;
                dim[n1]=ha_set[n].size;
                begadd[n1]=ha_set[n].begadd;
                break;
              }
            }*/
            if (strcmp(varset[n1],ha_set[ha_cof[i].setid[n1]].setname)==0) {
              dims=dims*ha_set[ha_cof[i].setid[n1]].size;
              dim[n1]=ha_set[ha_cof[i].setid[n1]].size;
              begadd[n1]=ha_set[ha_cof[i].setid[n1]].begadd;
              supsetid[n1]=0;
            } else {
              for (n=0; n<nset; n++) {
                if (strcmp(varset[n1],ha_set[n].setname)==0) {
                  dims=dims*ha_set[n].size;
                  dim[n1]=ha_set[n].size;
                  begadd[n1]=ha_set[n].begadd;
                  for(l1=0; l1<MAXSUPSET; l1++)if(strcmp(ha_set[ha_set[n].subsetid[l1]].setname,ha_set[ha_cof[i].setid[n1]].setname)==0) {
                      supsetid[n1]=l1;
                      break;
                    }
                  break;
                }
              }
            }
          }
          //printf("dim %d\n",dims);
          /*if (ha_cof[i].size>1) {
            dim[ha_cof[i].size-1]=dim1[1];
            varsize[ha_cof[i].size-1]=ha_cof[i].dims[1];
            begadd[ha_cof[i].size-1]=begadd1[1];
            //printf("dim %d\n",dim[ha_var[i].size-1]);
            dim[ha_cof[i].size-2]=dim1[0];
            varsize[ha_cof[i].size-2]=ha_cof[i].dims[0];
            begadd[ha_cof[i].size-2]=begadd1[0];
            //printf("dim %d\n",dim[ha_cof[i].size-2]);
            n=1;
            for (n1=ha_cof[i].size-3; n1>-1; n1--) {
              dim[n1]=dim1[ha_cof[i].size-n];
              varsize[n1]=ha_cof[i].dims[ha_cof[i].size-n];
              begadd[n1]=begadd1[ha_cof[i].size-n];
              n++;
              //printf("dim %d\n",dim[ha_cof[i].size-n]);
            }
          }*/
          /*dim[0]=dim1[1];
          varsize[0]=ha_cof[i].dims[1];
          begadd[0]=begadd1[1];
          dim[1]=dim1[0];
          varsize[1]=ha_cof[i].dims[0];
          begadd[1]=begadd1[0];
          for (n1=2; n1<ha_cof[i].size; n1++) {
          dim[n1]=dim1[n1];
          varsize[n1]=ha_cof[i].dims[n1];
          begadd[n1]=begadd1[n1];
          }*/

          antidim[0]=dim[1];
          //antidim1[0]=ha_cof[i].dims[1];
          antidim[1]=1;
          //antidim1[1]=1;
          antidim[2]=antidim[0]*dim[0];
          //antidim1[2]=antidim1[0]*ha_cof[i].dims[0];
          //printf("anti %d dim %d\n",antidim[0],dim[0]);
          //printf("anti %d dim %d\n",antidim[1],dim[1]);
          //printf("anti %d dim %d\n",antidim[2],dim[2]);
          if (ha_cof[i].size>3) {
            for (n1=3; n1<ha_cof[i].size; n1++) {
              antidim[n1]=antidim[n1-1]*dim[n1-1];
              //antidim1[n1]=antidim1[n1-1]*ha_cof[i].dims[n1-1];
              //printf("anti %d dim %d\n",antidim[n1],dim[n1]);
            }
          }
          //strcpy(vname,readitem);
          filehandle1 = fopen(iodata[k0].filname,"r");
          if(filehandle1==NULL){
            printf("Error!!! No such %s file!\n",iodata[k0].filname);
            return -1;
          }
          while (fgets(line,DATREADLINE,filehandle1)) {
            readitem = strtok(line,"\"");
            readitem = strtok(NULL,"\"");
            //printf("read %s\n",readitem);
            if (readitem != NULL) {
              strcpy(vname1,readitem);              //printf("read %s\n",readitem);
              n1=0;
              while (vname1[n1]!='\0'){
                if(vname1[n1]==' '){
                  vname1[n1]='\0';
                  break;
                }
                n1++;
              }
              //printf("vname %s\n",header);
              if (strncmp(vname1,header,strlen(header)) == 0&&strlen(header)==strlen(vname1)) {
                printf("dim %s\n",readitem);
                count2=2;
                recount=0;
                while (fgets(line,DATREADLINE,filehandle1)) {
                  count1=0;
                  copyvar=line;
                  while (*copyvar==' ' || *copyvar=='\t') {
                    copyvar++;  // skip whitespaces
                  }
                  if (*copyvar=='\r') {
                    copyvar++;
                  }
                  if (*copyvar=='\n') {
                    count1=1;
                  }
                  if (count1!=1) {
                    n=ha_cgenchf(line,',');
                    if (n==0) {
                      readitem = strtok(line,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if (ha_cof[i].size>1) {
                        for (n1=ha_cof[i].size-1; n1>1; n1--) {
                          index[n1]=(uvadd) l1/antidim[n1];
                          //printf("indx %d\n",index[n1]);
                          l1=l1-antidim[n1]*index[n1];
                        }
                        index[0]=(uvadd) l1/antidim[0];
                        l1=l1-antidim[0]*index[0];
                        index[1]=(uvadd) l1/antidim[1];
                        l1=l1-antidim[1]*index[1];
                      } else {
                        index[0]=recount;
                      }
                      recount1=0;
                      for (n1=0; n1<ha_cof[i].size; n1++) {
                        recount1=recount1+ha_setele[begadd[n1]+index[n1]].setsh[supsetid[n1]]*ha_cof[i].antidims[n1];
                        //printf("reacoun %d ssh %d\n",recount1,ha_setele[begadd[n1]+index[n1]].setsh);
                      }
                      //printf("reacoun %d ssh %d\n",recount1,ha_setele[recount1].setsh);
                      ha_cofele[ha_cof[i].begadd+recount1].cofval=val;
                      recount++;
                      if (recount>=dims) {
                        break;
                      }
                    } else {
                      for (n2=0; n2<n; n2++) {
                        //printf("recount %d\n",recount);
                        if (n2==0) {
                          readitem = strtok(line,",");
                        } else {
                          readitem = strtok(NULL,",");
                        }
                        val=atof(readitem);//sscanf(readitem, "%lf", &val);
                        l1=recount;
                        if (ha_cof[i].size>1) {
                          for (n1=ha_cof[i].size-1; n1>1; n1--) {
                            index[n1]=(uvadd) l1/antidim[n1];
                            //printf("indx %d\n",index[n1]);
                            l1=l1-antidim[n1]*index[n1];
                          }
                          index[0]=(uvadd) l1/antidim[0];
                          l1=l1-antidim[0]*index[0];
                          index[1]=(uvadd) l1/antidim[1];
                          l1=l1-antidim[1]*index[1];
                          //printf("indx %d\n",index[1]);
                          //printf("indx %d\n",index[0]);
                          //printf("recount %d\n",recount);
                        } else {
                          index[0]=recount;
                        }
                        recount1=0;
                        for (n1=0; n1<ha_cof[i].size; n1++) {
                          recount1=recount1+ha_setele[begadd[n1]+index[n1]].setsh[supsetid[n1]]*ha_cof[i].antidims[n1];
                          //printf("reacoun %d ssh %d\n",recount1,ha_setele[begadd[n1]+index[n1]].setsh);
                        }
                        //printf("reacoun1 %d\n",recount1);
                        ha_cofele[ha_cof[i].begadd+recount1].cofval=val;
                        //newrecount=rdatconv(ha_cof,i,recount);
                        //printf("newrec %d rec %d\n",newrecount,recount);
                        recount++;
                      }
                      readitem = strtok(NULL,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if (ha_cof[i].size>1) {
                        for (n1=ha_cof[i].size-1; n1>1; n1--) {
                          index[n1]=(uvadd) l1/antidim[n1];
                          //printf("indx %d\n",index[n1]);
                          l1=l1-antidim[n1]*index[n1];
                        }
                        index[0]=(uvadd) l1/antidim[0];
                        l1=l1-antidim[0]*index[0];
                        index[1]=(uvadd) l1/antidim[1];
                        l1=l1-antidim[1]*index[1];
                      } else {
                        index[0]=recount;
                      }
                      recount1=0;
                      for (n1=0; n1<ha_cof[i].size; n1++) {
                        recount1=recount1+ha_setele[begadd[n1]+index[n1]].setsh[supsetid[n1]]*ha_cof[i].antidims[n1];
                        //printf("reacoun %d ssh %d\n",recount1,ha_setele[begadd[n1]+index[n1]].setsh);
                      }
                      //printf("reacoun %d\n",recount1);
                      ha_cofele[ha_cof[i].begadd+recount1].cofval=val;
                      //newrecount=rdatconv(ha_cof,i,recount);
                      //printf("newrec %d rec %d\n",newrecount,recount);
                      recount++;
                      if (recount>=dims) {
                        break;
                      }
                    }
                  }

                }
              }

            }
            if(feof(filehandle1)) {
              break;
            }
          }
          if(count2!=2){
            printf("Error!!!! %s is not a header in data file!\n",header);
            return -1;
          }
          fclose(filehandle1);
          break;
        }
      }
      if(i<ncof)count2=1;
      for (i=0; i<nvar; i++) {
        vname1= strtok(ha_var[i].cofname,"(");
        //printf("vname %s varname %s varsize %d\n",vname,vname1,ha_var[i].size);
        if (strcmp(vname1,vname)==0) { //,vsize ha_var[i].cofname
          dims=1;
          for (n1=0; n1<ha_var[i].size; n1++) {
            /*for (n=0; n<nset; n++) {
              if (strcmp(varset[n1],ha_set[n].setname)==0) {
                dims=dims*ha_set[n].size;
                dim[n1]=ha_set[n].size;
                begadd[n1]=ha_set[n].begadd;
                break;
              }
            }*/
            if (strcmp(varset[n1],ha_set[ha_var[i].setid[n1]].setname)==0) {
              dims=dims*ha_set[ha_var[i].setid[n1]].size;
              dim[n1]=ha_set[ha_var[i].setid[n1]].size;
              begadd[n1]=ha_set[ha_var[i].setid[n1]].begadd;
              supsetid[n1]=0;
            } else {
              for (n=0; n<nset; n++) {
                if (strcmp(varset[n1],ha_set[n].setname)==0) {
                  dims=dims*ha_set[n].size;
                  dim[n1]=ha_set[n].size;
                  begadd[n1]=ha_set[n].begadd;
                  for(l1=0; l1<MAXSUPSET; l1++)if(strcmp(ha_set[ha_set[n].subsetid[l1]].setname,ha_set[ha_var[i].setid[n1]].setname)==0) {
                      supsetid[n1]=l1;
                      break;
                    }
                  break;
                }
              }
            }
          }
          //printf("dim %d\n",dims);
          /*if (ha_var[i].size>1) {
            dim[ha_var[i].size-1]=dim1[1];
            varsize[ha_var[i].size-1]=ha_var[i].dims[1];
            begadd[ha_var[i].size-1]=begadd1[1];
            //printf("dim %d\n",dim[ha_var[i].size-1]);
            dim[ha_var[i].size-2]=dim1[0];
            varsize[ha_var[i].size-2]=ha_var[i].dims[0];
            begadd[ha_var[i].size-2]=begadd1[0];
            //printf("dim %d\n",dim[ha_var[i].size-2]);
            n=1;
            for (n1=ha_var[i].size-3; n1>-1; n1--) {
              dim[n1]=dim1[ha_var[i].size-n];
              varsize[n1]=ha_var[i].dims[ha_var[i].size-n];
              begadd[n1]=begadd1[ha_var[i].size-n];
              n++;
              //printf("dim %d\n",dim[ha_var[i].size-n]);
            }
          }*/
          /*dim[0]=dim1[1];
          varsize[0]=ha_var[i].dims[1];
          begadd[0]=begadd1[1];
          dim[1]=dim1[0];
          varsize[1]=ha_var[i].dims[0];
          begadd[1]=begadd1[0];
          for (n1=2; n1<ha_var[i].size; n1++) {
          dim[n1]=dim1[n1];
          varsize[n1]=ha_var[i].dims[n1];
          begadd[n1]=begadd1[n1];
          }*/

          antidim[0]=dim[1];
          //antidim1[0]=ha_var[i].dims[1];
          antidim[1]=1;
          //antidim1[1]=1;
          antidim[2]=antidim[0]*dim[0];
          //antidim1[2]=antidim1[0]*ha_var[i].dims[0];
          //printf("anti %d dim %d\n",antidim[0],dim[0]);
          //printf("anti %d dim %d\n",antidim[1],dim[1]);
          //printf("anti %d dim %d\n",antidim[2],dim[2]);
          if (ha_var[i].size>3) {
            for (n1=3; n1<ha_var[i].size; n1++) {
              antidim[n1]=antidim[n1-1]*dim[n1-1];
              //antidim1[n1]=antidim1[n1-1]*ha_var[i].dims[n1-1];
              //printf("anti %d dim %d\n",antidim[n1],dim[n1]);
            }
          }
          //strcpy(vname,readitem);
          filehandle1 = fopen(iodata[k0].filname,"r");
          if(filehandle1==NULL){
            printf("Error!!! No such %s file!\n",iodata[k0].filname);
            return -1;
          }
          while (fgets(line,DATREADLINE,filehandle1)) {
            readitem = strtok(line,"\"");
            readitem = strtok(NULL,"\"");
            //printf("read %s\n",readitem);
            if (readitem != NULL) {
              strcpy(vname1,readitem);              //printf("read %s\n",readitem);
              n1=0;
              while (vname1[n1]!='\0'){
                if(vname1[n1]==' '){
                  vname1[n1]='\0';
                  break;
                }
                n1++;
              }
              //printf("read %s\n",readitem);
              //printf("vname %s\n",header);
              if (strncmp(vname1,header,strlen(header)) == 0&&strlen(header)==strlen(vname1)) {
                //printf("dim %s\n",readitem);
                count2=4;
                recount=0;
                while (fgets(line,DATREADLINE,filehandle1)) {
                  count1=0;
                  copyvar=line;
                  while (*copyvar==' ' || *copyvar=='\t') {
                    copyvar++;  // skip whitespaces
                  }
                  if (*copyvar=='\r') {
                    copyvar++;
                  }
                  if (*copyvar=='\n') {
                    count1=1;
                  }
                  if (count1!=1) {
                    n=ha_cgenchf(line,',');
                    if (n==0) {
                      readitem = strtok(line,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if (ha_var[i].size>1) {
                        for (n1=ha_var[i].size-1; n1>1; n1--) {
                          index[n1]=(uvadd) l1/antidim[n1];
                          //printf("indx %d\n",index[n1]);
                          l1=l1-antidim[n1]*index[n1];
                        }
                        index[0]=(uvadd) l1/antidim[0];
                        l1=l1-antidim[0]*index[0];
                        index[1]=(uvadd) l1/antidim[1];
                        l1=l1-antidim[1]*index[1];
                      } else {
                        index[0]=recount;
                      }
                      recount1=0;
                      for (n1=0; n1<ha_var[i].size; n1++) {
                        recount1=recount1+ha_setele[begadd[n1]+index[n1]].setsh[supsetid[n1]]*ha_var[i].antidims[n1];
                        //printf("reacoun %d ssh %d\n",recount1,ha_setele[begadd[n1]+index[n1]].setsh);
                      }
                      //printf("reacoun %d ssh %d\n",recount1,ha_setele[recount1].setsh);
                      ha_varele[ha_var[i].begadd+recount1].cofval=val;
                      recount++;
                      if (recount>=dims) {
                        break;
                      }
                    } else {
                      for (n2=0; n2<n; n2++) {
                        //printf("recount %d\n",recount);
                        if (n2==0) {
                          readitem = strtok(line,",");
                        } else {
                          readitem = strtok(NULL,",");
                        }
                        val=atof(readitem);//sscanf(readitem, "%lf", &val);
                        l1=recount;
                        if (ha_var[i].size>1) {
                          for (n1=ha_var[i].size-1; n1>1; n1--) {
                            index[n1]=(uvadd) l1/antidim[n1];
                            //printf("indx %d\n",index[n1]);
                            l1=l1-antidim[n1]*index[n1];
                          }
                          index[0]=(uvadd) l1/antidim[0];
                          l1=l1-antidim[0]*index[0];
                          index[1]=(uvadd) l1/antidim[1];
                          l1=l1-antidim[1]*index[1];
                          //printf("indx %d\n",index[1]);
                          //printf("indx %d\n",index[0]);
                          //printf("recount %d\n",recount);
                        } else {
                          index[0]=recount;
                        }
                        recount1=0;
                        for (n1=0; n1<ha_var[i].size; n1++) {
                          recount1=recount1+ha_setele[begadd[n1]+index[n1]].setsh[supsetid[n1]]*ha_var[i].antidims[n1];
                          //printf("reacoun %d ssh %d\n",recount1,ha_setele[begadd[n1]+index[n1]].setsh);
                        }
                        //printf("reacoun1 %d\n",recount1);
                        ha_varele[ha_var[i].begadd+recount1].cofval=val;
                        //newrecount=rdatconv(ha_var,i,recount);
                        //printf("newrec %d rec %d\n",newrecount,recount);
                        recount++;
                      }
                      readitem = strtok(NULL,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if (ha_var[i].size>1) {
                        for (n1=ha_var[i].size-1; n1>1; n1--) {
                          index[n1]=(uvadd) l1/antidim[n1];
                          //printf("indx %d\n",index[n1]);
                          l1=l1-antidim[n1]*index[n1];
                        }
                        index[0]=(uvadd) l1/antidim[0];
                        l1=l1-antidim[0]*index[0];
                        index[1]=(uvadd) l1/antidim[1];
                        l1=l1-antidim[1]*index[1];
                      } else {
                        index[0]=recount;
                      }
                      recount1=0;
                      for (n1=0; n1<ha_var[i].size; n1++) {
                        recount1=recount1+ha_setele[begadd[n1]+index[n1]].setsh[supsetid[n1]]*ha_var[i].antidims[n1];
                        //printf("reacoun %d ssh %d\n",recount1,ha_setele[begadd[n1]+index[n1]].setsh);
                      }
                      //printf("reacoun %d\n",recount1);
                      ha_varele[ha_var[i].begadd+recount1].cofval=val;
                      //newrecount=rdatconv(ha_var,i,recount);
                      //printf("newrec %d rec %d\n",newrecount,recount);
                      recount++;
                      if (recount>=dims) {
                        break;
                      }
                    }
                  }

                }
              }

            }
            if(feof(filehandle1)) {
              break;
            }
          }
          if(count2!=4){
            printf("Error!!!! %s is not a header in data file!\n",header);
            return -1;
          }
          fclose(filehandle1);
          break;
        }
      }
      if(count2==0&&i==ncof){
        printf("Error!!!! %s is not a variable coefficient or parameter!\n",vname);
        return -1;
      }
    }
  }
  fclose(filehandle);
  return j;
}

uvadd rdatconv(hcge_cof *ha_var,uvadd varindx,uvadd loops,ha_cgeset *ha_set) {
  uvadd dcount,dcount1[MAXVARDIM],bdcount1[MAXVARDIM],cdcount1[MAXVARDIM],l;//i,
  if(ha_var[varindx].size<2) {
    printf("Variable dimension should bigger than 1");
    return -1;
  }
  for (dcount=0; dcount<ha_var[varindx].size; dcount++) {
    dcount1[dcount]=(long int) loops/ha_set[ha_var[varindx].setid[dcount]].size;
    loops=loops-dcount1[dcount]*ha_set[ha_var[varindx].setid[dcount]].size;
    //printf("dc %d dims %d loops %d\n",dcount1[dcount],ha_var[varindx].dims[dcount],loops);
  }
  for (dcount=1; dcount<ha_var[varindx].size; dcount++) {
    l=ha_var[varindx].size-dcount;
    bdcount1[l]=dcount1[dcount];
    cdcount1[l-2]=ha_set[ha_var[varindx].setid[l-1]].size;
  }
  bdcount1[0]=dcount1[ha_var[varindx].size-2];
  bdcount1[1]=dcount1[ha_var[varindx].size-1];
  cdcount1[ha_var[varindx].size-1]=1;
  cdcount1[0]=ha_set[ha_var[varindx].setid[ha_var[varindx].size-1]].size;
  l=0;
  for (dcount=0; dcount<ha_var[varindx].size; dcount++) {
    l=l+bdcount1[dcount]*cdcount1[dcount];
  }
  return l;
}

uvadd hcge_dsum(char *formulain, char *commsyntax, hcge_sumcof *sum_cof,ha_cgesetindx *arSet,ha_cgeset *ha_set,uvdim nset,uvdim fdim,int j) {
  char *readitem,*p,*p1,*p2,interchar2[TABREADLINE],argu[TABREADLINE];//,line5[TABREADLINE]
  char interchar[TABREADLINE],interchar1[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],line3[TABREADLINE],line4[TABREADLINE],tempname[NAMESIZE];
  uvadd i=0,k=0,k1=0,length,ncur=0,ncuri,l,l1,l2,l3,l4,l5,l6,l7,sup;
  length=strlen(formulain);
  readitem=formulain;
  //printf("readitem %s\n",readitem);
  //printf("readitem %s\n",commsyntax);
  while (i<length) {
    k=ha_cgefind(readitem,commsyntax);
    //printf("k %d l %d\n",k,length);
    //printf("readitem %s\n",readitem);
    if (k==-1) {
      //printf("sumcount1 %d\n",j);
      return 0;
    }
    if (k==0) {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      strcpy(line1,readitem);
      ha_cgecutsum(line);
      k1=ha_cgefind(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      } else {
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
        //strcpy(ha_set[sum_cof[j].sumsetid].setname,p);
        //printf("sunset %s\n",p);
        for (l7=0; l7<nset; l7++) if(strcmp(p,ha_set[l7].setname)==0) {
            //sum_cof[j].sumsize=ha_set[l7].size;
            //sum_cof[j].ssetbegadd=ha_set[l7].begadd;
            sum_cof[j].sumsetid=l7;
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
                    //strcpy(sum_cof[j].dimsets[l3],arSet[l5].arSet);
                    //sum_cof[j].dims[l3]=arSet[l5].SetSize;
                    sum_cof[j].setid[l3]=arSet[l5].setid;
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
                  //strcpy(ha_set[sum_cof[j].setid[l3]].setname,p1);
                  for (l7=0; l7<nset; l7++) if(strcmp(p1,ha_set[l7].setname)==0) {
                      //sum_cof[j].dims[l3]=ha_set[l7].size;
                      sum_cof[j].setid[l3]=l7;
                      //for(sup=0;sup<MAXSUPSET;sup++)sum_cof[j].subsetid[l3][sup]=ha_set[l7].subsetid[sup];
                      //sum_cof[j].supsetsize[l3]=ha_set[l7].supersetsize;
                      //sum_cof[j].dimssetbegadd[l3]=ha_set[l7].begadd;
                      break;
                    }
                }
                l3++;
              }
            }
          } else {
            for (l1=0; l1<l; l1++) {
              if(l1==0) {
                p = strtok(argu,",");
              } else {
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
                      //strcpy(sum_cof[j].dimsets[l3],arSet[l5].arSet);
                      //printf("set %s\n",arSet[l5].arSet);
                      //sum_cof[j].dims[l3]=arSet[l5].SetSize;
                      sum_cof[j].setid[l3]=arSet[l5].setid;
                      //for(sup=0;sup<MAXSUPSET;sup++)sum_cof[j].subsetid[l3][sup]=arSet[l5].subsetid[sup];
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
                    //strncpy(sum_cof[j].dimsets[l3],p1,p2-p1);
                    //sum_cof[j].dimsets[l3][p2-p1]='\0';
                    strncpy(tempname,p1,p2-p1);
                    tempname[p2-p1]='\0';
                    //printf("setname %s\n",sum_cof[j].dimsets[l3]);
                    for (l7=0; l7<nset; l7++) if(strcmp(tempname,ha_set[l7].setname)==0) {
                        //sum_cof[j].dims[l3]=ha_set[l7].size;
                        sum_cof[j].setid[l3]=l7;
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
        } else {
          if (interchar[strlen(interchar)-1]=='{') {
            interchar[strlen(interchar)-1]='\0';
          } else {
            strcat(interchar,"}");
          }
        }
        while(ha_cgefrstr(formulain,line4,interchar));
        //printf("readitem %s\n",formulain);
        //printf("line4 %s inter %s\n",line4,interchar);
        sum_cof[j].size=l3;
        return 1;
      }
    } else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      strcpy(line1,readitem);
      ha_cgecutsum(line);
      //printf("readitem %s\n",line);
      k1=ha_cgefind(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      } else {
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
        //strcpy(tempname,p);
        //printf("sunset %s\n",p);
        for (l7=0; l7<nset; l7++) if(strcmp(p,ha_set[l7].setname)==0) {
            //sum_cof[j].sumsize=ha_set[l7].size;
            //sum_cof[j].ssetbegadd=ha_set[l7].begadd;
            sum_cof[j].sumsetid=l7;
            //printf("sumsetname %s\n",ha_set[sum_cof[j].sumsetid].setname);
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
                strcat(interchar,p);
                l6=0;
                for (l5=0; l5<fdim-1; l5++) if(strcmp(p,arSet[l5].arIndx)==0) {
                    //strcpy(sum_cof[j].dimsets[l3],arSet[l5].arSet);
                    //sum_cof[j].dims[l3]=arSet[l5].SetSize;
                    sum_cof[j].setid[l3]=arSet[l5].setid;
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
                  //strcpy(tempname,p1);
                  for (l7=0; l7<nset; l7++) if(strcmp(p1,ha_set[l7].setname)==0) {
                      //sum_cof[j].dims[l3]=ha_set[l7].size;
                      sum_cof[j].setid[l3]=l7;
                      //sum_cof[j].subsetid[l3]=ha_set[l7].subsetid;
                      //sum_cof[j].supsetsize[l3]=ha_set[l7].supersetsize;
                      //sum_cof[j].dimssetbegadd[l3]=ha_set[l7].begadd;
                      break;
                    }
                }
                l3++;
              }
            }
          } else {
            for (l1=0; l1<l; l1++) {
              if(l1==0) {
                p = strtok(argu,",");
              } else {
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
                  //printf("p %s sumindx %s j %d\n",p,sum_cof[j].sumindx,j);
                  strcat(interchar,sum_cof[j].dimnames[l3]);
                  strcat(interchar,",");
                  l6=0;
                  for (l5=0; l5<fdim-1; l5++) if(strcmp(p,arSet[l5].arIndx)==0) {
                      //strcpy(sum_cof[j].dimsets[l3],arSet[l5].arSet);
                      //sum_cof[j].dims[l3]=arSet[l5].SetSize;
                      sum_cof[j].setid[l3]=arSet[l5].setid;
                      //printf("p %s sumindx %s j %d l3 %d\n",p,ha_set[sum_cof[j].setid[l3]].setname,j,l3);
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
                    //strncpy(sum_cof[j].dimsets[l3],p1,p2-p1);
                    //sum_cof[j].dimsets[l3][p2-p1]='\0';
                    strncpy(tempname,p1,p2-p1);
                    tempname[p2-p1]='\0';
                    for (l7=0; l7<nset; l7++) if(strcmp(tempname,ha_set[l7].setname)==0) {
                        //sum_cof[j].dims[l3]=ha_set[l7].size;
                        sum_cof[j].setid[l3]=l7;
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
        } else {
          if (interchar[strlen(interchar)-1]=='{') {
            interchar[strlen(interchar)-1]='\0';
          } else {
            strcat(interchar,"}");
          }
        }
        while(ha_cgefrstr(formulain,line4,interchar));
        //printf("readitem %s\n",formulain);
        //printf("line4 %s inter %s\n",line4,interchar);
        sum_cof[j].size=l3;
        //for (l5=0; l5<sum_cof[j].size; l5++)printf("4p %s sumindx %s j %d l3 %d\n",p,ha_set[sum_cof[j].setid[l5]].setname,j,l5);
        return 1;
      }
    } else {
      i=i+k+4;
      readitem=formulain+i;
    }
  }
  return 0;
}
/*
int hcge_rsum(char *formulain, char *commsyntax, hcge_sumcof *sum_cof)
{
  char *readitem,*p;
  char interchar[NAMESIZE],interchar1[NAMESIZE],line[TABREADLINE];
  uvadd i=0,k=0,length;
  int j=0;
  length=strlen(formulain);
  readitem=formulain;
  //printf("readitem %s\n",readitem);
  //printf("readitem %s\n",commsyntax);
  while (i<length) {
    k=ha_cgefind(readitem,commsyntax);
    //printf("k %d l %d\n",k,length);
    //printf("readitem %s\n",readitem);
    if (k==-1) {
      break;
    }
    if (k==0) {
      i=i+k+4;
      readitem=formulain+i;
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
      if (j>10000) {
        strcat(interchar,"ha_cgesum");
        printf("Error: Too many sum\n");
      }
      strcat(interchar,interchar1);
      strcpy(sum_cof[j].sumname,interchar);
      strcpy(line,readitem);
      p = strtok(line,",");
      strcpy(sum_cof[j].sumindx,p);
      p = strtok(NULL,",");
      strcpy(sum_cof[j].sumset,p);
      j++;
    } else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      i=i+k+4;
      readitem=formulain+i;
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
      if (j>10000) {
        strcat(interchar,"ha_cgesum");
        printf("Error: Too many sum\n");
      }
      strcat(interchar,interchar1);
      strcpy(sum_cof[j].sumname,interchar);
      strcpy(line,readitem);
      p = strtok(line,",");
      strcpy(sum_cof[j].sumindx,p);
      p = strtok(NULL,",");
      strcpy(sum_cof[j].sumset,p);
      j++;
    } else {
      i=i+k+4;
      readitem=formulain+i;
    }
  }
  return j;
}*/

int hcge_nsum(char *formulain, char *commsyntax) {
  char *readitem;
  int j=0,i=0,k=0,length;
  length=strlen(formulain);
  readitem=formulain;
  //printf("readitem %s l %d\n",readitem,length);
  //printf("readitem %s\n",commsyntax);
  while (i<length) {
    k=ha_cgefind(readitem,commsyntax);
    //printf("k %d l %d\n",k,length);
    if (k==-1) {
      break;
    }
    if (k==0) {
      j++;
      i=i+k+4;
      readitem=formulain+i;
    } else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==','||formulain[i+k-1]=='='||formulain[i+k-1]=='['||formulain[i+k-1]=='{'||formulain[i+k-1]==' ') {
      j++;
      i=i+k+4;
      readitem=formulain+i;
      //printf("read %s\n",readitem);
    } else {
      i=i+k+4;
      readitem=formulain+i;
      //printf("read %s\n",readitem);
    }
  }
  return j;
}

uvadd ha_cgencof(char *fname, char *commsyntax) {
  FILE * filehandle;
  char line[TABREADLINE]="\0";
  //char *commsyntax="coefficient";
  uvadd j=0;
  filehandle = fopen(fname,"r");
  while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
    //printf("line %s\n",line);
    if (strstr(line,"(default")==NULL) {
      j++;
    }
  }
  fclose(filehandle);
  return j;
}

uvadd hcge_rexo(char *fname, char *commsyntax,ha_cgeexovar *ha_cgeshock, hcge_cof *ha_var,uvadd nvar,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele) {
  FILE * filehandle;
  char line[TABREADLINE]="\0",*readitem=NULL,*p=NULL,*p1=NULL,vname[TABREADLINE],argu[TABREADLINE];//,linecopy[TABREADLINE]
  uvadd i,j,n=0,l=0,l1,l2,n1,m,dims,dcount,exodims,supsetid[MAXSUPSET],sup,doublepr[MAXSUPSET];
  int k1,k2;
  bool check;
  filehandle = fopen(fname,"r");
  //hcge_cof *exo_var= (hcge_cof *) calloc (1,sizeof(hcge_cof));
  while (ha_cgercls(commsyntax,filehandle,line)) {
    //printf("read %s\n",line);
    //while (ha_cgefrstr(line,"\n", ""));
    while (ha_cgedrcmt(line,"!"));
    ha_cgefrstr(line,";", " ;");
    //while (ha_cgefrstr(line,"  ", " "));
    //while (ha_cgefrstr(line,"\n", " "));
    while (ha_cgefrstr(line,"\n", " "));
    while (ha_cgefrstr(line,"\r", " "));
    while (ha_cgefrstr1(line," p_", " "));
    while (ha_cgefrstr1(line," c_", " "));
    while (ha_cgefrstr(line,"  ", " "));
    while (ha_cgefrstr(line," ,", ","));
    while (ha_cgefrstr(line,", ", ","));
    while (ha_cgefrstr(line," )", ")"));
    while (ha_cgefrstr(line," (", "("));
    while (ha_cgefrstr(line," \"", "\""));
    while (ha_cgefrstr(line,"\" ", "\""));
    //while (ha_cgefrstr(line,"\n", ""));
    //readitem = strtok(line," ");
    k1=0;
    k2=0;
    //printf("read %s\n",line);
    while (line[k1]!= '\0') {
      if(line[k1]=='\"') {
        if(k2==0) {
          k2=1;
        } else {
          k2=0;
        }
      } else {
        if(k2==0) {
          line[k1]=tolower((int)line[k1]);
        }
      }
      k1++;
    }
    n1=ha_cgenchf(line,' ');
    readitem = strchr(line,' ');
    readitem++;
    //printf("read %s\n n1 %d\n",readitem,n1);
    //while (ha_cgefrstr(line,"\n", " "));
    //printf("read %s\n n1 %d\n",readitem,n1);
    for (i=0; i<n1-1; i++) {
      p = strchr(readitem,' ');
      strncpy(vname,readitem,p-readitem);
      vname[p-readitem]='\0';
      //printf("vname %s\n",vname);
      readitem=p;
      readitem++;
      //printf("read %s\n",readitem);
      check=true;
      if (strchr(vname,'(')==NULL) {
        for (j=0; j<nvar; j++) if (strcmp(vname,ha_var[j].cofname)==0) {
            dims=ha_var[j].matsize;
            //printf("var %s dims %d n %d\n",ha_var[j].cofname,dims,n);
            //for (m=0; m<ha_var[j].size; m++) {
            //  dims=dims*ha_var[j].dims[m];
            //}
            if (dims==0) {
              ha_cgeshock[ha_var[j].begadd].ShockId=true;
              n=n+1;
            } else for (l=0; l<dims; l++) {
                ha_cgeshock[ha_var[j].begadd+l].ShockId=true;
                n=n+1;
              }
            check=false;
            //printf("n %d\n",n);
            break;
          }
          if(j==nvar&&vname[0]!=';')printf("Error!!! There is no such var %s\n",vname);
      } else {
        dims=1;
        p = strtok(vname,"(");
        //printf("var %s n %d\n",vname,n);
        for(sup=0; sup<MAXSUPSET; sup++)supsetid[sup]=0;
        for (j=0; j<nvar; j++) if (strcmp(p,ha_var[j].cofname)==0) {
            ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (ha_var[j].size,sizeof(ha_cgesetindx));
            uvadd *exoantidim= (uvadd *) calloc (ha_var[j].size,sizeof(uvadd));
            switch (ha_var[j].size) {
            case 1:
              p = strtok(NULL,")");//p = strtok(vname,")");
              if (strchr(p,'"')!=NULL) {
                strcpy(argu,p);
                p=strchr(argu,'"');
                p+=1;
                p1=strchr(p,'"');
                *p1='\0';
                //p = strtok(argu,"\"");
                //p = strtok(NULL,"\"");
                //printf("KHKNKJ %s\n",p);
                for (l1=0; l1<ha_set[ha_var[j].setid[0]].size; l1++)
                  if (strcmp(p,ha_setele[ha_set[ha_var[j].setid[0]].begadd+l1].setele)==0) {
                    ha_cgeshock[ha_var[j].begadd+l1].ShockId=true;
                    n=n+1;
                    check=false;
                    break;
                  }
                  if(l1==ha_set[ha_var[j].setid[0]].size)printf("Error!!! There is no element %s in the set %s in cof/var %s\n",p,ha_set[ha_var[j].setid[0]].setname,ha_var[j].cofname);
              } else {
                if(strcmp(p,ha_set[ha_var[j].setid[0]].setname)==0) {
                  //arSet[0].SetSize=ha_var[j].dims[0];
                  //arSet[0].SetBegAdd=ha_var[j].setbegadd1[0];
                  arSet[0].setid=ha_var[j].setid[0];
                  supsetid[0]=0;
                  exoantidim[0]=1;//arSet[0].SetSize;
                  dims=ha_set[arSet[0].setid].size;
                } else {
                  for (l1=0; l1<nset; l1++) if (strcmp(p,ha_set[l1].setname)==0) {
                      //arSet[0].SetSize=ha_set[l1].size;
                      //arSet[0].SetBegAdd=ha_set[l1].begadd;
                      //printf("set name %s ",p);
                      //for(sup=0; sup<MAXSUPSET; sup++)printf("sup %i subid %i\n",sup,ha_set[l1].subsetid[sup]);
                      arSet[0].setid=l1;
                      for(sup=1; sup<MAXSUPSET; sup++)if(ha_set[l1].subsetid[sup]==ha_var[j].setid[0])supsetid[0]=sup;
                      exoantidim[0]=1;//arSet[0].SetSize;
                      dims=ha_set[arSet[0].setid].size;
                      //printf("subid %d dims %d l1 %d\n",supsetid[0],dims,l1);
                      break;
                    }
                    if(l1==nset)printf("Error!!! There is no set %s in cof/var %s\n",p,ha_var[j].cofname);
                }
              }
              break;
            default:
              p = strtok(NULL,")");
              strcat(p,",");
              //printf("p %s\n",p);
              for (l=0; l<ha_var[j].size; l++) {
                if (l==0) {
                  p = strtok(p,",");
                } else {
                  p = strtok(NULL,",");
                }
                //printf("p %s\n",p);
                if (*p=='"') {
                  strcpy(argu,p+1);
                  ha_cgefrchr(argu,'"','\0');
                  for(l1=0; l1<ha_set[ha_var[j].setid[l]].size; l1++) if (strcmp(argu,ha_setele[ha_set[ha_var[j].setid[l]].begadd+l1].setele)==0) {
                      //printf("argu %s\n",argu);
                      //arSet[l].SetSize=1;
                      //arSet[l].SetBegAdd=ha_var[j].setbegadd1[l]+l1;//Special case
                      //arSet[l].setid=
                      doublepr[l]=l1;
                      supsetid[l]=-2;
                      break;
                    }
                    if(l1==ha_set[ha_var[j].setid[l]].size)printf("Error!!! There is no element %s in the set %s in cof/var %s\n",argu,ha_set[ha_var[j].setid[l]].setname,ha_var[j].cofname);
                  //dims=dims;
                } else {
                  if(strcmp(p,ha_set[ha_var[j].setid[l]].setname)==0) {
                    arSet[l].setid=ha_var[j].setid[l];
                    //supsetid[l]=0;
                    //for(sup=0;sup<MAXSUPSET;sup++)if(ha_set[l1].subsetid[sup]==l1)supsetid[0]=sup;
                    dims=dims*ha_set[arSet[l].setid].size;
                  } else {
                    for(l1=0; l1<nset; l1++) if (strcmp(p,ha_set[l1].setname)==0) {
                        //arSet[l].SetSize=ha_set[l1].size;
                        //printf("set %d\n",arSet[l].SetSize);
                        //arSet[l].SetBegAdd=ha_set[l1].begadd;
                        //arSet[l].subsetid=ha_set[l1].subsetid;
                        arSet[l].setid=l1;
                        for(sup=1; sup<MAXSUPSET; sup++)if(ha_set[l1].subsetid[sup]==ha_var[j].setid[l])supsetid[l]=sup;
                        dims=dims*ha_set[arSet[l].setid].size;
                        break;
                      }
                      if(l1==nset)printf("Error!!! There is no set %s in cof/var %s\n",p,ha_var[j].cofname);
                  }
                }
                //printf("dims %d\n",dims);
              }
              exoantidim[ha_var[j].size-1]=1;
              for (l1=ha_var[j].size-2; l1>-1; l1--) {
                if(supsetid[l1+1]>-2)exoantidim[l1]=exoantidim[l1+1]*ha_set[arSet[l1+1].setid].size;
                else exoantidim[l1]=exoantidim[l1+1];
              }
              break;
            }
            //for (l1=0;l1<ha_var[j].size;l1++) printf("exoantidims %d dims %d\n",exoantidim[l1],dims);
            if(check) for (l=0; l<dims; l++) {
                m=l;
                exodims=0;
                for (dcount=0; dcount<ha_var[j].size; dcount++) {
                  l1=(uvadd) m/exoantidim[dcount];
                  m=m-l1*exoantidim[dcount];
                  if(supsetid[dcount]==-2) {
                    l2=doublepr[dcount];//ha_set[arSet[dcount].setid].begadd;  //Special case
                  } else {
                    l2=ha_setele[ha_set[arSet[dcount].setid].begadd+l1].setsh[supsetid[dcount]];
                  }
                  //printf("l1 %d l2 %d set beg %d\n",l1,l2,ha_set[arSet[dcount].setid].begadd);
                  exodims=exodims+l2*ha_var[j].antidims[dcount];
                }
                //printf("ha_var %s bedg %d exo %d\n",ha_var[j].cofname,ha_var[j].begadd,exodims);
                ha_cgeshock[ha_var[j].begadd+exodims].ShockId=true;
                n=n+1;
              }
            free(arSet);
            free(exoantidim);
            break;
          }
          if(j==nvar)printf("Error!!! There is no such var %s\n",p);
        //printf("n %d\n",n);
      }
      //readitem = strtok(NULL," ");
    }
  }
  //free(exo_var);
  fclose(filehandle);
  return n;
}

uvadd hcge_rshock(char *fname, char *commsyntax,ha_cgeexovar *ha_cgeshock,uvadd nvarele, hcge_cof *ha_var,uvadd nvar,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele,uvdim subints) {
  FILE * filehandle;
  char line[DATREADLINE],linecopy[DATREADLINE],argu[TABREADLINE];
  line[0]='\0';
  linecopy[0]='\0';
  argu[0]='\0';
  char *readitem=NULL,*p=NULL;
  uvdim varnset;
  uvdim k1,k2,dimindx[MAXVARDIM],dimbegadd[MAXVARDIM],dimssize[MAXVARDIM],arSet[MAXVARDIM],antidim[MAXVARDIM];
  uvadd j,l=0,dims,n1,l1,l2,dcount,supsetid[MAXSUPSET],sup;
  ha_cgetype val;
      if ( (filehandle = fopen(fname,"r")) == NULL ) {
        printf("Error opening file %s\n",fname);
      }
  //filehandle = fopen(fname,"r");
  //printf("OK!!!\n");
  while (ha_cgertabl(commsyntax,filehandle,line,DATREADLINE)) {
    //printf("OK!!!\n");
    //printf("read %s\n",line);
    //printf("OK!!!!!\n");
    //while (ha_cgefrstrvbz(line,"\n", " ",DATREADLINE));
    ha_chrfrall(line,'\r',' ');
    ha_chrfrall(line,'\n',' ');
    //printf("OK!!!!!\n");
    //while (ha_cgedrcmt(line,"!"));
    //printf("read1\n %s\n",line);
    k1=0;
    k2=0;
    //printf("read1 %s\n",line);
    //printf("OK!!!!!\n");
    while (line[k1]!= '\0') {
      //if(line[k1]=='\"') {
      //  if(k2==0) {
      //    k2=1;
      //  } else {
      //    k2=0;
      //  }
      //} else {
      //  if(k2==0) {
      line[k1]=tolower((int)line[k1]);
      //  }
      //}
      k1++;
    }
    //printf("k1 %d\n",k1);
    if(k1>=DATREADLINE){
      printf("Error!!!! The shock file line (from Shock to ;) is too long!\n");
      return -1;
    }
    //printf("read2 %s\n",line);
    ha_cgefrstrvbz(line,";", " ;",DATREADLINE);
    while (ha_cgefrstr(line,"= ", "="));
    while (ha_cgefrstrvbz(line,"  ", " ",DATREADLINE));
    while (ha_cgefrstrvbz1(line," p_", " ",DATREADLINE));
    while (ha_cgefrstrvbz1(line," c_", " ",DATREADLINE));
    //printf("shock %s\n",line);
    varnset=ha_cgenchf(line,',');
    j=ha_cgenchf(line,'(');
    for(sup=0; sup<MAXSUPSET; sup++)supsetid[sup]=-1;
    if (j>0) {
      varnset=varnset+1;
    }
    //printf("line %s\n",line);
    if(varnset==0){
      while(ha_cgefrstr(line," =", "="));
    }

    if (line[0]==' ') {
      readitem = strtok(line+1," ");
    } else {
      readitem = strtok(line," ");
    }
    //printf("line %s\n",readitem);
    if (varnset==0) {
      readitem = strtok(NULL,"=");
      //printf("readitem %s\n",readitem);
      for (j=0; j<nvar; j++) if (strcmp(readitem,ha_var[j].cofname)==0) {
          //readitem = strtok(NULL," ");
          readitem = strtok(NULL,";");
          //printf("readitem %s\n",readitem);
          k1=ha_cgerevfind(readitem,"uniform");
          //printf("readitem %s\n",readitem+k1+1);
          if(k1!=-1){
            ha_cgeshock[ha_var[j].begadd].ShockVal=atof(readitem+k1+1)/subints;
          }else{
            ha_cgeshock[ha_var[j].begadd].ShockVal=atof(readitem)/subints;
          }
          l=l+1;
          break;
        }
    } else {
      readitem = strtok(NULL,"(");
      //printf("read1 %s\n",readitem);
      //printf("line %s\n",line);
      for (j=0; j<nvar; j++) {
        if (strcmp(readitem,ha_var[j].cofname)==0) {
          dims=1;
          readitem = strtok(NULL,")");
          strcpy(argu,readitem);
          readitem = strtok(NULL,"=");
          //readitem = strtok(NULL," ");
          readitem = strtok(NULL,";");
          //printf("read %s\n",readitem);
          strcpy(linecopy,readitem);
          strcat(linecopy," ");
          strcat(argu,",");
          while (ha_cgefrstr(argu," ", ""));
          for (n1=0; n1<MAXVARDIM; n1++) {
            dimindx[n1]=-1;
            dimbegadd[n1]=-1;
            dimssize[n1]=0;
            arSet[n1]=0;
            antidim[n1]=0;
          }
          //printf("OK!!!\n");
          for (n1=0; n1<ha_var[j].size; n1++) {
            if(n1==0) {
              p=strtok(argu,",");
            } else {
              p=strtok(NULL,",");
            }
            if(strchr(p,'\"')==NULL) {
              //printf("set %s\n",p);
              if((strcmp(ha_set[ha_var[j].setid[n1]].setname,p)==0)) {
                //printf("p %s var %s\n",p,ha_set[ha_var[j].setid[n1]].setname);
                dimbegadd[n1]=ha_set[ha_var[j].setid[n1]].begadd;//ha_set[k1].begadd;
                dimssize[n1]=ha_set[ha_var[j].setid[n1]].size;//ha_set[k1].size;
                dims=dims*dimssize[n1];
                supsetid[n1]=0;
              } else {
                for (k1=0; k1<nset; k1++)if (strcmp(ha_set[k1].setname,p)==0) {
                    dimbegadd[n1]=ha_set[k1].begadd;
                    dimssize[n1]=ha_set[k1].size;
                    dims=dims*ha_set[k1].size;
                    for(sup=0; sup<MAXSUPSET; sup++)if(ha_set[k1].subsetid[sup]==ha_var[j].setid[n1]) {
                        supsetid[n1]=sup;
                        break;
                      }
                    if(sup==MAXSUPSET){
                      printf("Error!!!! %s in variable %s in shock file is not a correct set at that position!\n",p,ha_var[j].cofname);
                      return -1;
                    }
                    break;
                  }
                  if(k1==nset){
                    printf("Error!!!! %s in variable %s in shock file is not a set!\n",p,ha_var[j].cofname);
                    return -1;
                  }
                //printf("p %s var %s\n",p,ha_set[ha_var[j].setid[n1]].setname);
              }
            } else {
              k2=0;
              while (p[k2] != '\0') {
                k2++;
              }
              for (k1=0; k1<ha_set[ha_var[j].setid[n1]].size; k1++)if (strncmp(ha_setele[ha_set[ha_var[j].setid[n1]].begadd+k1].setele,p+1,k2-2)==0&&ha_setele[ha_set[ha_var[j].setid[n1]].begadd+k1].setele[k2-2]=='\0') {
                  dimindx[n1]=ha_setele[ha_set[ha_var[j].setid[n1]].begadd+k1].setsh[0];
                  dimssize[n1]=1;
                  supsetid[n1]=-2;
                  //printf("dimindx %d p %s k2 %d\n",dimindx[n1],p,k2);
                  break;
                }
              if(k1==ha_set[ha_var[j].setid[n1]].size){
                printf("Error!!!! %s in variable %s in shock file is not a set element!\n",p,ha_var[j].cofname);
                return -1;
              }
              //printf("p %s var %s\n",p,ha_set[ha_var[j].setid[n1]].setname);
            }
          }
          //printf("OK!!!\n");
          if (ha_var[j].size>0)antidim[0]=dimssize[1];
          if (ha_var[j].size>1)antidim[1]=1;
          if (ha_var[j].size>2)antidim[2]=antidim[0]*dimssize[0];
          if (ha_var[j].size>3) {
            for (n1=3; n1<ha_var[j].size; n1++) {
              antidim[n1]=antidim[n1-1]*dimssize[n1-1];
            }
          }
          /*for (n1=ha_var[j].size; n1>0; n1--) {
            if(n1==ha_var[j].size) {
              antidim[n1-1]=1;
            } else {
              antidim[n1-1]=antidim[n1]*dimssize[n1];
            }
            //printf("dim %d anti %d\n",dimssize[n1-1],antidim[n1-1]);
          }*/
          //strcpy(linecopy,readitem);
          k1=ha_cgerevfind(linecopy,"uniform");
          //printf("OK1!!!\n");
          if(k1>-1) {
            val=atof(linecopy+k1+1);//sscanf(readitem, "%lf", &val);
            //printf("ha_var %s dims %d vval %lf\n",ha_var[j].cofname,dims,val);
            for (n1=0; n1<dims; n1++) {
              l2=n1;
              if(ha_var[j].size>1) {
                for (dcount=ha_var[j].size-1; dcount>1; dcount--) {
                  if(dimindx[dcount]>-1) {
                    //l1=dimindx[dcount];
                    arSet[dcount]=dimindx[dcount];//ha_setele[dimbegadd[dcount]+l1].setsh;
                  } else {
                    l1=(uvadd) l2/antidim[dcount];
                    if(supsetid[dcount]>-2)arSet[dcount]=ha_setele[dimbegadd[dcount]+l1].setsh[supsetid[dcount]];
                    else arSet[dcount]=dimindx[dcount];
                    l2=l2-l1*antidim[dcount];
                  }
                }
                if(dimindx[0]>-1) {
                  arSet[0]=dimindx[0];
                } else {
                  l1=(uvadd) l2/antidim[0];
                  if(supsetid[0]>-2)arSet[0]=ha_setele[dimbegadd[0]+l1].setsh[supsetid[0]];
                  else arSet[0]=ha_setele[dimbegadd[0]+l1].setsh[0];//dimindx[0];
                  l2=l2-l1*antidim[0];
                }
                if(dimindx[1]>-1) {
                  arSet[1]=dimindx[1];
                } else {
                  l1=(uvadd) l2/antidim[1];
                  if(supsetid[1]>-2)arSet[1]=ha_setele[dimbegadd[1]+l1].setsh[supsetid[1]];
                  else arSet[1]=ha_setele[dimbegadd[1]+l1].setsh[0];//l1;//dimindx[1];
                  //arSet[1]=ha_setele[dimbegadd[1]+l1].setsh[supsetid[dcount]];
                }
              } else {
                arSet[0]=l2;
              }
              l2=0;
              for (dcount=0; dcount<ha_var[j].size; dcount++) {
                l2=l2+arSet[dcount]*ha_var[j].antidims[dcount];
              }
              ha_cgeshock[ha_var[j].begadd+l2].ShockVal=val/subints;
              l=l+1;
            }
            //for (n1=0; n1<dims; n1++) {
              //ha_cgeshock[ha_var[j].begadd+n1].ShockVal=val;
              //l=l+1;
            //}
          } else {
            //printf("line %s\n",linecopy);
            for (n1=0; n1<dims; n1++) {
              //printf("n1 %d dims %d\n",n1,dims);
              if(n1==0) {
                if(linecopy[0]==' ') {
                  readitem=strtok(linecopy+1," ");
                  //printf("line %s\n",linecopy);
                  //readitem=strtok(NULL," ");
                } else {
                  readitem=strtok(linecopy," ");
                }
              } else {
                readitem=strtok(NULL," ");
              }
              val=atof(readitem);
              //printf("n1 %ld read %s\n",n1,readitem);
              l2=n1;
              if(ha_var[j].size>1) {
                for (dcount=ha_var[j].size-1; dcount>1; dcount--) {
                  if(dimindx[dcount]>-1) {
                    //l1=dimindx[dcount];
                    arSet[dcount]=dimindx[dcount];//ha_setele[dimbegadd[dcount]+l1].setsh;
                  } else {
                    l1=(uvadd) l2/antidim[dcount];
                    if(supsetid[dcount]>-2)arSet[dcount]=ha_setele[dimbegadd[dcount]+l1].setsh[supsetid[dcount]];
                    else arSet[dcount]=dimindx[dcount];
                    l2=l2-l1*antidim[dcount];
                  }
                  //index[dcount]=(uvadd) l1/antidim[dcount];
                  //l1=l1-antidim[dcount]*index[dcount];
                }
                if(ha_var[j].size==2){//compliance
                if(dimindx[1]>-1) {
                  arSet[1]=dimindx[1];
                } else {
                  l1=(uvadd) l2/dimssize[0];
                  //printf("ds %d\n",dimssize[0]);
                  if(supsetid[1]>-2)arSet[1]=ha_setele[dimbegadd[1]+l1].setsh[supsetid[1]];
                  else arSet[1]=ha_setele[dimbegadd[1]+l1].setsh[0];//dimindx[1];
                  l2=l2-l1*dimssize[0];
                }
                if(dimindx[0]>-1) {
                  arSet[0]=dimindx[0];
                } else {
                  l1=(uvadd) l2;
                  if(supsetid[0]>-2)arSet[0]=ha_setele[dimbegadd[0]+l1].setsh[supsetid[0]];
                  else arSet[0]=ha_setele[dimbegadd[0]+l1].setsh[0];//l1;//dimindx[0];
                  //arSet[0]=ha_setele[dimbegadd[0]+l1].setsh[supsetid[0]];
                }
                }else{
                if(dimindx[0]>-1) {
                  //l1=dimindx[dcount];
                  arSet[0]=dimindx[0];//ha_setele[dimbegadd[dcount]+l1].setsh;
                } else {
                  l1=(uvadd) l2/antidim[0];
                  if(supsetid[0]>-2)arSet[0]=ha_setele[dimbegadd[0]+l1].setsh[supsetid[0]];
                  else arSet[0]=ha_setele[dimbegadd[0]+l1].setsh[0];//dimindx[0];
                  l2=l2-l1*antidim[0];
                }
                if(dimindx[1]>-1) {
                  arSet[1]=dimindx[1];//ha_setele[dimbegadd[dcount]+l1].setsh;
                } else {
                  l1=(uvadd) l2/antidim[1];
                  if(supsetid[1]>-2)arSet[1]=ha_setele[dimbegadd[1]+l1].setsh[supsetid[1]];
                  else arSet[1]=ha_setele[dimbegadd[1]+l1].setsh[0];//dimindx[0];
                  //arSet[1]=ha_setele[dimbegadd[1]+l1].setsh[supsetid[dcount]];
                  //l2=l2-l1*antidim[1];
                }
                //index[0]=(uvadd) l1/antidim[0];
                //l1=l1-antidim[0]*index[0];
                //printf("reacoun1 %d\n",recount1);
                //index[1]=(uvadd) l1/antidim[1];
                //l1=l1-antidim[1]*index[1];
                }
              } else {
                if(supsetid[0]>-2)arSet[0]=ha_setele[dimbegadd[0]+l2].setsh[supsetid[0]];
                else arSet[0]=dimindx[0];
                //printf("aaaaaaaaaaaa %d\n",arSet[0]);
                //arSet[0]=l2;
              }
              l2=0;
              for (dcount=0; dcount<ha_var[j].size; dcount++) {
                l2=l2+arSet[dcount]*ha_var[j].antidims[dcount];
                //printf("reacoun %d ssh %d\n",recount1,ha_setele[begadd[n1]+index[n1]].setsh);
              }
              //printf("aaaaaaaaaaaa %d var %s begadd %d exo %d\n",l2,ha_var[j].cofname,ha_var[j].begadd,ha_var[j].begadd+l2);
              //printf("reacoun1 %d\n",recount1);
              //ha_varele[ha_var[j].begadd+l2].cofval=val;
              /*for (dcount=0; dcount<ha_var[j].size; dcount++) {
                //printf("dcount %d\n",dcount);
                if(dimindx[dcount]>-1) {
                  //l1=dimindx[dcount];
                  arSet[dcount]=dimindx[dcount];//ha_setele[dimbegadd[dcount]+l1].setsh;
                } else {
                  l1=(uvadd) l2/antidim[dcount];
                  arSet[dcount]=ha_setele[dimbegadd[dcount]+l1].setsh;
                  l2=l2-l1*antidim[dcount];
                }
              }
              l2=0;
              //printf("n1 %d read %s\n",n1,readitem);
              for (dcount=0; dcount<ha_var[j].size; dcount++) {
                l2=l2+arSet[dcount]*ha_var[j].antidims[dcount];
                //if (strcmp("a1tot",ha_var[j].cofname)==0&&dcount==1)printf("i %d arset %d\n",dcount,arSet[dcount]);
              }*/
              ha_cgeshock[ha_var[j].begadd+l2].ShockVal=val/subints;
              l=l+1;
            }
          }
          break;
        }
      }
      if(j==nvar){
        printf("Error!!!! %s in shock file is not a variable!\n",readitem);
        return -1;
      }
    }
  }
  fclose(filehandle);
  return l;
}

int hcge_defvar(char *fname, hcge_cof *record, uvadd ncof) {
  char line[TABREADLINE];
  FILE * filehandle;
  char commsyntax[]="variable";
  filehandle = fopen(fname,"r");
  uvadd i;
  while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
    if(ha_cgefind(line,"default")>0) {
      while (ha_cgefrstr(line," ", ""));
      if(ha_cgefind(line,"default=levels")>0) for (i=0; i<ncof; i++) {
          record[i].level_par=true;
        }
      if(ha_cgefind(line,"default=change")>0) for (i=0; i<ncof; i++) {
          record[i].change_real=true;
        }
    }
  }
  fclose(filehandle);
  return 1;
}

uvadd hcge_rvar(char *fname, char *commsyntax, hcge_cof *record, uvadd ncof, ha_cgeset *ha_set,uvdim nset) {
  FILE * filehandle;
  char line[TABREADLINE]="\0",linecopy[TABREADLINE],setname1[NAMESIZE],setname[NAMESIZE],setname2[NAMESIZE],setname3[NAMESIZE],finditem[NAMESIZE],finditem1[NAMESIZE],finditem2[NAMESIZE],finditem3[NAMESIZE],vname[NAMESIZE];//,vnamecopy[NAMESIZE];
  char *tpnt=NULL;
  //char tempcofname[ncof][TABREADLINE];
  uvadd n,m,l,ncommsyntax=0,i,j=0,addi=0,add=0,orig;
  uvdim dcount;
  char *readitem=NULL;
  bool IsLevel=false,IsChange=false;
  while (commsyntax[ncommsyntax] != '\0') {
    ncommsyntax++;
  }
  filehandle = fopen(fname,"r");

  while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
    if(strstr(line,"(default=levels)")!=NULL) {
      IsLevel=true;
    }
    if(strstr(line,"(default=change)")!=NULL) {
      IsChange=true;
    }
    if(IsLevel) {
      record[j].level_par=true;
    }
    if(IsChange) {
      record[j].change_real=true;
    }
    /*while (ha_cgefrstr(line,"\n", ""));
    while (ha_cgefrstr(line,"\r", ""));
    while (ha_cgefrstr(line,"  ", " "));
    strcat(line, "\n");
    while (ha_cgedrcmt(line,"!"));
    while (ha_cgedrcmt(line,"#"));
    while (ha_cgefrstr(line,"  ", " "));*/
    readitem=strstr(line, "(ge ");
    //printf("line %s\n",line);
    if(readitem!=NULL){
      record[j].gltype=1;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      tpnt+=1;
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, "(ge ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
      //printf("line %s\n type %d val %lf read %s\n",linecopy,record[j].gltype,record[j].glval,readitem);
      //while (ha_cgefrstr(line," ", ""));
      //while (ha_cgefrstr(line,"(ge0)", ""));
      //printf("line %s\n type %d val %lf\n",line,record[j].gltype,record[j].glval);
    }
    readitem=strstr(line, "(gt ");
    if(readitem!=NULL){
      record[j].gltype=2;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      tpnt+=1;
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, "(gt ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, "(le ");
    if(readitem!=NULL){
      record[j].gltype=3;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      tpnt+=1;
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, "(le ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, "(lt ");
    if(readitem!=NULL){
      record[j].gltype=4;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      tpnt+=1;
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, "(lt ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }

    readitem=strstr(line, ",ge ");
    if(readitem!=NULL){
      record[j].gltype=1;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, ",ge ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
      printf("line %s\n type %d val %lf read %s\n",linecopy,record[j].gltype,record[j].glval,readitem);
      //while (ha_cgefrstr(line," ", ""));
      //while (ha_cgefrstr(line,"(ge0)", ""));
      printf("line %s\n type %d val %lf\n",line,record[j].gltype,record[j].glval);
    }
    readitem=strstr(line, ",gt ");
    if(readitem!=NULL){
      record[j].gltype=2;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, ",gt ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, ",le ");
    if(readitem!=NULL){
      record[j].gltype=3;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, ",le ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, ",lt ");
    if(readitem!=NULL){
      record[j].gltype=4;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, ",lt ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    while (ha_cgefrstr(line," ", ""));
    //while (ha_cgefrstr(line,"(ge0)", ""));
    //ha_cgefrstr1(line,"(PARAMETER)", "");
    if(ha_cgefind(line,"(change)")>0) {
      record[j].change_real=true;
      ha_cgefrstr1(line,"(change)", "");
    }
    if(ha_cgefind(line,"(percent_change)")>0) {
      record[j].change_real=false;
      ha_cgefrstr1(line,"(percent_change)", "");
    }
    if(ha_cgefind(line,"(levels)")>0) {
      record[j].level_par=true;
      ha_cgefrstr1(line,"(levels)", "");
    }
    if(ha_cgefind(line,"(linear)")>0) {
      record[j].level_par=false;
      ha_cgefrstr1(line,"(linear)", "");
    }
    if(ha_cgefind(line,"(linear,change)")>0) {
      record[j].change_real=true;
      record[j].level_par=false;
      ha_cgefrstr1(line,"(linear,change)", "");
    }
    if(ha_cgefind(line,"(change,linear)")>0) {
      record[j].change_real=true;
      record[j].level_par=false;
      ha_cgefrstr1(line,"(change,linear)", "");
    }
    strcpy(linecopy,line);
    //printf("%s\n",line);
    if (strstr(line,"(default")==NULL) {
      n=ha_cgenchf(line,')');
      orig=0;
      if (strstr(line,"(orig_level")!=NULL) orig=1;
      //printf("n %d\n",n);
      if (n-orig>1) {
        for (i=1; i<n-orig; i++) {
          if (i==1) {
            readitem = strtok(line,",");
          } else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
          finditem[0]='(';
          finditem[1]='\0';
          finditem1[0]=',';
          finditem1[1]='\0';
          finditem2[0]=',';
          finditem2[1]='\0';
          strcat(finditem,readitem);
          strcat(finditem1,readitem);
          strcat(finditem2,readitem);
          strcpy(finditem3,finditem);
          //strcpy(finditem,readitem);
          //strcpy(finditem1,readitem);
          //strcpy(finditem2,readitem);
          //strcat("(\0",finditem);
          strcat(finditem,",");
          //strcat(",\0",finditem1);
          strcat(finditem1,",");
          //strcat(",\0",finditem2);
          strcat(finditem2,")");
          strcat(finditem3,")");
          //printf("%s\n",finditem);
          //printf("%s\n",finditem1);
          //printf("%s\n",finditem2);
          readitem = strtok(NULL,")");
          setname[0]='(';
          setname[1]='\0';
          setname1[0]=',';
          setname1[1]='\0';
          setname2[0]=',';
          setname2[1]='\0';
          strcat(setname,readitem);
          strcat(setname1,readitem);
          strcat(setname2,readitem);
          strcpy(setname3,setname);
          //strcpy(setname,readitem);
          //strcpy(setname1,readitem);
          //strcpy(setname2,readitem);
          //strcat("(\0",setname);
          strcat(setname,",");
          //strcat(",\0",setname1);
          strcat(setname1,",");
          //strcat(",\0",setname2);
          strcat(setname2,")");
          strcat(setname3,")");
          while (ha_cgefrstr(linecopy, finditem, setname));
          while (ha_cgefrstr(linecopy, finditem1, setname1));
          while (ha_cgefrstr(linecopy, finditem2, setname2));
          while (ha_cgefrstr(linecopy, finditem3, setname3));
        }
        //printf("linecopy %s\n",linecopy);
        for (i=1; i<n; i++) {
          if (i==1) {
            readitem = strtok(linecopy,")");
          } else {
            readitem = strtok(NULL,")");
          }
        }
        //printf("vname %s\n",readitem);
        readitem = strtok(NULL,";");
        //printf("snem %s j %d ncof %d\n",readitem,j,ncof);
        readitem = strtok(readitem,"(");
        strcpy(record[j].cofname,readitem);
        readitem = strtok(NULL,")");
        strcpy(vname,readitem);
        strcat(vname,",");
        dcount=ha_cgenchf(vname,',');
        add=1;
        //printf("vname1 %s\n",vname);
        for (m=0; m<dcount; m++) {
          if(m==0) {
            readitem = strtok(vname,",");
          } else {
            readitem = strtok(NULL,",");
          }
          for (l=0; l<nset; l++) {
            //printf("setname %s\t vname %s\n",ha_set[l].setname,vname);
            if (strcmp(ha_set[l].setname,readitem)==0) {
              //strcpy(record[j].dimsets[m],readitem);
              //record[j].dims[m]=ha_set[l].size;
              add=add*ha_set[l].size;
              record[j].setid[m]=l;
              //printf("j %d l %d setname %s\t vname %s setid \n",j,l,ha_set[l].setname,record[j].cofname);
              //printf("size %d\n",record[i].size);
              break;
            }
          }
        }
        record[j].begadd=addi;
        record[j].size=dcount;
        record[j].matsize=add;
        record[j].antidims[dcount-1]=1;
        for (l=record[j].size-2; l>-1; l--) {
          record[j].antidims[l]=record[j].antidims[l+1]*ha_set[record[j].setid[l+1]].size;
        }
        addi=addi+add;
        //printf("addi %d add %d\n",addi,add);
        //strcpy(tempcofname[j].line,readitem);
        //for (i=0;i<strlen(readitem);i++)tempcofname[j][i]=readitem[i];
        //tempcofname[j][strlen(readitem)]='\0';
        //tempcofname[j]);
      } else {
        if (n-orig==1) {
          printf("%s\n","Syntax error! Not enough parences");
        } else {
          if(orig==1) {
            strcpy(setname,line+ncommsyntax);
            ha_cgefrstr(setname,";", "");
            ha_cgefrstr(setname,"\n", "");
            readitem=strchr(setname,')');
            readitem++;
            strcpy(record[j].cofname,readitem);
            //printf("snem %s j %d ncof %d\n",record[j].cofname,j,ncof);
            record[j].begadd=addi;
            record[j].size=0;
            record[j].matsize=1;
            addi=addi+1;
          } else {
            strcpy(setname,line+ncommsyntax);
            ha_cgefrstr(setname,";", "");
            ha_cgefrstr(setname,"\n", "");
            //printf("snem %s j %d ncof %d\n",setname,j,ncof);
            strcpy(record[j].cofname,setname);
            record[j].begadd=addi;
            record[j].size=0;
            record[j].matsize=1;
            addi=addi+1;
          }
          //strcpy(tempcofname[j].line,setname);
        }
      }
      j++;
    }
  }
  fclose(filehandle);
  /*
    j=0;
    for (i=0; i<ncof; i++) {
      strcpy(vname,tempcofname[i]);
      //strcpy(vnamecopy,vname);
      dcount=ha_cgenchf(vname,')');
      //printf("vname %s\n",vname);
      if (dcount==0) {
        //strcpy(vname,vnamecopy);
        strcpy(record[i].cofname,vname);
        record[i].begadd=j;
        record[i].size=0;
        record[i].matsize=1;
        j=j+1;
      }
      if (dcount==1) {
        //strcpy(vname,vnamecopy);
        dcount=ha_cgenchf(vname,',');
        if (dcount==0) {
          readitem = strtok(vname,"(");
          strcpy(record[i].cofname,readitem);
          readitem = strtok(NULL,")");
          strcpy(vname,readitem);
          m=0;
          while (vname[m] != '\0') {
            m++;
          }
          for (l=0; l<nset; l++) {
            //printf("setname %s\t vname %s\n",ha_set[l].setname,vname);
            if (strncmp(ha_set[l].setname,vname,m)==0) {
              record[i].begadd=j;
              record[i].size=1;
              record[i].dims[0]=ha_set[l].size;
              record[i].antidims[0]=1;
              record[i].matsize=ha_set[l].size;
              strcpy(record[i].dimsets[0],vname);
              j=j+ha_set[l].size;
              //printf("size %d\n",record[i].size);
              break;
            }
          }
        }
        if (dcount==1) {
          //printf("vname %s\n",vname);
          readitem = strtok(vname,"(");
          strcpy(record[i].cofname,readitem);
          readitem = strtok(NULL,",");
          strcpy(vname,readitem);
          m=0;
          while (vname[m] != '\0') {
            m++;
          }
          for (l=0; l<nset; l++) {
            //printf("setname %s\t vname %s\n",ha_set[l].setname,vname);
            if (strncmp(ha_set[l].setname,vname,m)==0) {
              strcpy(record[i].dimsets[0],vname);
              record[i].dims[0]=ha_set[l].size;
              //printf("size %d\n",record[i].size);
              break;
            }
          }
          readitem = strtok(NULL,")");
          strcpy(vname,readitem);
          m=0;
          while (vname[m] != '\0') {
            m++;
          }
          for (l=0; l<nset; l++) {
            //printf("setname %s\t vname %s\n",ha_set[l].setname,vname);
            if (strncmp(ha_set[l].setname,vname,m)==0) {
              strcpy(record[i].dimsets[1],vname);
              record[i].dims[1]=ha_set[l].size;
              //printf("size %d\n",record[i].size);
              break;
            }
          }
          record[i].begadd=j;
          record[i].size=2;
          record[i].antidims[1]=1;
          record[i].antidims[0]=record[i].dims[1];
          record[i].matsize=record[i].dims[0]*record[i].dims[1];
          j=j+record[i].matsize;
        }
        if (dcount>1) {
          //printf("vname %s\n",vname);
          add=1;
          readitem = strtok(vname,"(");
          strcpy(record[i].cofname,readitem);
          for(k=0; k<dcount; k++) {
            readitem = strtok(NULL,",");
            strcpy(vname,readitem);
            m=0;
            while (vname[m] != '\0') {
              m++;
            }
            for (l=0; l<nset; l++) {
              //printf("setname %s\t vname %s\n",ha_set[l].setname,vname);
              if (strncmp(ha_set[l].setname,vname,m)==0) {
                strcpy(record[i].dimsets[k],vname);
                record[i].dims[k]=ha_set[l].size;
                add=add*ha_set[l].size;
                //printf("size %d\n",record[i].size);
                break;
              }
            }
          }
          readitem = strtok(NULL,")");
          strcpy(vname,readitem);
          m=0;
          while (vname[m] != '\0') {
            m++;
          }
          for (l=0; l<nset; l++) {
            //printf("setname %s\t vname %s\n",ha_set[l].setname,vname);
            if (strncmp(ha_set[l].setname,vname,m)==0) {
              strcpy(record[i].dimsets[k],vname);
              record[i].dims[k]=ha_set[l].size;
              add=add*ha_set[l].size;
              //printf("size %d\n",record[i].size);
              break;
            }
          }
          record[i].begadd=j;
          record[i].size=dcount+1;
          record[i].matsize=add;
          record[i].antidims[dcount]=1;
          for (l=record[i].size-2; l>-1; l--) {
            record[i].antidims[l]=record[i].antidims[l+1]*record[i].dims[l+1];
          }
          j=j+add;
        }
      }
    }*/
  //for (i=0;i<ncof;i++) printf("oday %s\n",record[i].cofname);
  return addi;
}

uvadd ha_cgeralltime(ha_cgeset *ha_set,uvdim nset) {
  uvdim i;
  for(i=0; i<nset; i++) {
    if(ha_set[i].intertemp&&ha_set[i].subsetid[1]==-1) {
      break;
    }
  }
  if(i==nset) {
    return -1;
  } else {
    return i;
  }
}

uvadd hcge_req(char *fname, char *commsyntax, hcge_cof *record, uvadd ncof, ha_cgeset *ha_set,uvdim nset) {
  FILE * filehandle;
  char line[TABREADLINE]="\0",linecopy[TABREADLINE];
  uvadd n,i,j=0,ncommsyntax,l;
  char *readitem=NULL;
  filehandle = fopen(fname,"r");

  while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
    if (strstr(line,"(default")==NULL) {
      n=ha_cgenfind(line,"(all,");
      //printf("n %d\n",n);
      record[j].size=n;
      if (n>0) {
        readitem = strtok(line," ");
        readitem = strtok(NULL,"(");
        strcpy(record[j].cofname,readitem);
        ncommsyntax=0;
        while (record[j].cofname[ncommsyntax] != '\0') {
          ncommsyntax++;
        }
        if(record[j].cofname[ncommsyntax-1]==' ') {
          record[j].cofname[ncommsyntax-1]='\0';
        }
        for (i=0; i<n; i++) {
          //if(i>0)readitem = strtok(NULL,"(");
          readitem = strtok(NULL,",");
          readitem = strtok(NULL,",");
          readitem = strtok(NULL,")");
          //printf("read %s j %d\n",readitem,j);
          for(l=0; l<nset; l++)if(strcmp(ha_set[l].setname,readitem)==0) {
              record[j].setid[i]=l;
              break;
            }
        }
      } else {
        readitem = strtok(line," ");
        readitem = strtok(NULL," ");
        strcpy(record[j].cofname,readitem);
        //record[j].size=0;
        record[j].matsize=1;
      }
      j++;
    }
  }
  fclose(filehandle);
  return j;
}

uvadd hcge_rcof(char *fname, char *commsyntax, hcge_cof *record, uvadd ncof, ha_cgeset *ha_set,uvadd nset) {
  FILE * filehandle;
  char line[TABREADLINE]="\0",linecopy[TABREADLINE],setname1[TABREADLINE],setname[TABREADLINE],setname2[TABREADLINE],setname3[TABREADLINE],finditem[TABREADLINE],finditem1[TABREADLINE],finditem2[TABREADLINE],finditem3[TABREADLINE],vname[TABREADLINE];//,vnamecopy[NAMESIZE];
  //ha_line* tempcofname= (ha_line *) calloc(ncof,sizeof(ha_line *));
  //char tempcofname[ncof][TABLINESIZE];
  uvadd n,m,l,ncommsyntax=0,i=0,j=0,addi=0,add=0;
  uvdim dcount;
  char *readitem=NULL,*tpnt=NULL;
  while (commsyntax[ncommsyntax] != '\0') {
    ncommsyntax++;
  }
  filehandle = fopen(fname,"r");

  while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
    /*while (ha_cgefrstr(line,"\n", ""));
    while (ha_cgefrstr(line,"\r", ""));
    while (ha_cgefrstr(line,"  ", " "));
    strcat(line, "\n");
    while (ha_cgedrcmt(line,"!"));
    while (ha_cgedrcmt(line,"#"));
    while (ha_cgefrstr(line,"  ", " "));*/
    readitem=strstr(line, "(ge ");
    //printf("line %s\n",line);
    if(readitem!=NULL){
      record[j].gltype=1;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      tpnt+=1;
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, "(ge ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
      //printf("line %s\n type %d val %lf read %s\n",linecopy,record[j].gltype,record[j].glval,readitem);
      //while (ha_cgefrstr(line," ", ""));
      //while (ha_cgefrstr(line,"(ge0)", ""));
      //printf("line %s\n type %d val %lf\n",line,record[j].gltype,record[j].glval);
    }
    readitem=strstr(line, "(gt ");
    if(readitem!=NULL){
      record[j].gltype=2;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      tpnt+=1;
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, "(gt ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, "(le ");
    if(readitem!=NULL){
      record[j].gltype=3;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      tpnt+=1;
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, "(le ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, "(lt ");
    if(readitem!=NULL){
      record[j].gltype=4;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      tpnt+=1;
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, "(lt ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }

    readitem=strstr(line, ",ge ");
    if(readitem!=NULL){
      record[j].gltype=1;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, ",ge ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, ",gt ");
    if(readitem!=NULL){
      record[j].gltype=2;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, ",gt ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, ",le ");
    if(readitem!=NULL){
      record[j].gltype=3;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, ",le ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, ",lt ");
    if(readitem!=NULL){
      record[j].gltype=4;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, ",lt ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    
    while (ha_cgefrstr(line," ", ""));
    //while (ha_cgefrstr(line,"(ge0)", ""));
    //while (ha_cgefrstr(line,"(ge 0)", ""));
    ha_cgefrstr1(line,"parameter", "");
    ha_cgefrstr1(line,"change", "");
    ha_cgefrstr1(line,"integer", "");
    ha_cgefrstr1(line,"()", "");
    ha_cgefrstr1(line,"(,)", "");
    strcpy(linecopy,line);
    //printf("%s\n",line);
    if (strstr(line,"(default")==NULL) {
      n=ha_cgenchf(line,')');
      //printf("n %d\n",n);
      if (n>1) {
        for (i=1; i<n; i++) {
          if (i==1) {
            readitem = strtok(line,",");
          } else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
          finditem[0]='(';
          finditem[1]='\0';
          finditem1[0]=',';
          finditem1[1]='\0';
          finditem2[0]=',';
          finditem2[1]='\0';
          strcat(finditem,readitem);
          strcat(finditem1,readitem);
          strcat(finditem2,readitem);
          strcpy(finditem3,finditem);
          //strcpy(finditem,readitem);
          //strcpy(finditem1,readitem);
          //strcpy(finditem2,readitem);
          //strcat("(\0",finditem);
          strcat(finditem,",");
          //strcat(",\0",finditem1);
          strcat(finditem1,",");
          //strcat(",\0",finditem2);
          strcat(finditem2,")");
          strcat(finditem3,")");
          //printf("%s\n",finditem);
          //printf("%s\n",finditem1);
          //printf("%s\n",finditem2);
          readitem = strtok(NULL,")");
          setname[0]='(';
          setname[1]='\0';
          setname1[0]=',';
          setname1[1]='\0';
          setname2[0]=',';
          setname2[1]='\0';
          strcat(setname,readitem);
          strcat(setname1,readitem);
          strcat(setname2,readitem);
          strcpy(setname3,setname);
          //strcpy(setname,readitem);
          //strcpy(setname1,readitem);
          //strcpy(setname2,readitem);
          //strcat("(\0",setname);
          strcat(setname,",");
          //strcat(",\0",setname1);
          strcat(setname1,",");
          //strcat(",\0",setname2);
          strcat(setname2,")");
          strcat(setname3,")");
          while (ha_cgefrstr(linecopy, finditem, setname));
          while (ha_cgefrstr(linecopy, finditem1, setname1));
          while (ha_cgefrstr(linecopy, finditem2, setname2));
          while (ha_cgefrstr(linecopy, finditem3, setname3));
        }
        //printf("linecopy %s\n",linecopy);
        for (i=1; i<n; i++) {
          if (i==1) {
            readitem = strtok(linecopy,")");
          } else {
            readitem = strtok(NULL,")");
          }
        }
        //printf("vname %s\n",readitem);
        readitem = strtok(NULL,";");
        //printf("snem %s j %d ncof %d\n",readitem,j,ncof);
        readitem = strtok(readitem,"(");
        strcpy(record[j].cofname,readitem);
        if(record[j].cofname[0]=='c'&&record[j].cofname[1]=='_'){
          printf("Errors!!!! c_ and C_ are reserved for change variable!!!\nPlease change the name of coefficient %s\n",record[j].cofname);
          return -1;
        }
        readitem = strtok(NULL,")");
        strcpy(vname,readitem);
        strcat(vname,",");
        dcount=ha_cgenchf(vname,',');
        add=1;
        for (m=0; m<dcount; m++) {
          if(m==0) {
            readitem = strtok(vname,",");
          } else {
            readitem = strtok(NULL,",");
          }
          for (l=0; l<nset; l++) {
            //printf("setname %s\t vname %s\n",ha_set[l].setname,vname);
            if (strcmp(ha_set[l].setname,readitem)==0) {
              //strcpy(record[j].dimsets[m],readitem);
              //record[j].dims[m]=ha_set[l].size;
              add=add*ha_set[l].size;
              record[j].setid[m]=l;
              //printf("size %d\n",record[i].size);
              break;
            }
          }
        }
        record[j].begadd=addi;
        record[j].size=dcount;
        record[j].matsize=add;
        record[j].antidims[dcount-1]=1;
        for (l=record[j].size-2; l>-1; l--) {
          record[j].antidims[l]=record[j].antidims[l+1]*ha_set[record[j].setid[l+1]].size;
        }
        addi=addi+add;
        //strcpy(tempcofname[j].line,readitem);
        //for (i=0;i<strlen(readitem);i++)tempcofname[j][i]=readitem[i];
        //tempcofname[j][strlen(readitem)]='\0';
        //tempcofname[j]);
      } else {
        if (n==1) {
          printf("%s\n","Syntax error! Not enough parences");
          printf("line %s\n",line);
        } else {
          strcpy(setname,line+ncommsyntax);
          ha_cgefrstr(setname,";", "");
          ha_cgefrstr(setname,"\n", "");
          //printf("snem %s j %d ncof %d\n",setname,j,ncof);
          strcpy(record[j].cofname,setname);
          record[j].begadd=addi;
          record[j].size=0;
          record[j].matsize=1;
          addi=addi+1;
          //printf("cofname %s j %d\n",setname,j);
          //strcpy(tempcofname[j].line,setname);
        }
      }
      j++;
    }
  }
  fclose(filehandle);

  //j=0;
  /*for (i=0; i<ncof; i++) {
    strcpy(vname,tempcofname[i].line);
    //strcpy(vnamecopy,vname);
    dcount=ha_cgenchf(vname,')');
    //printf("vname %s\n",vname);
    if (dcount==0) {
      //strcpy(vname,vnamecopy);
      strcpy(record[i].cofname,vname);
      record[i].begadd=j;
      record[i].size=0;
      record[i].matsize=1;
      j=j+1;
    }
    if (dcount==1) {
      //strcpy(vname,vnamecopy);
      dcount=ha_cgenchf(vname,',');
      if (dcount==0) {
        readitem = strtok(vname,"(");
        strcpy(record[i].cofname,readitem);
        readitem = strtok(NULL,")");
        strcpy(vname,readitem);
        m=0;
        while (vname[m] != '\0') {
          m++;
        }
        for (l=0; l<nset; l++) {
          //printf("setname %s\t vname %s\n",ha_set[l].setname,vname);
          if (strncmp(ha_set[l].setname,vname,m)==0) {
            record[i].begadd=j;
            record[i].size=1;
            record[i].dims[0]=ha_set[l].size;
            record[i].antidims[0]=1;
            record[i].matsize=ha_set[l].size;
            strcpy(record[i].dimsets[0],vname);
            j=j+ha_set[l].size;
            //printf("size %d\n",record[i].size);
            break;
          }
        }
      }
      if (dcount==1) {
        //printf("vname %s\n",vname);
        readitem = strtok(vname,"(");
        strcpy(record[i].cofname,readitem);
        readitem = strtok(NULL,",");
        strcpy(vname,readitem);
        m=0;
        while (vname[m] != '\0') {
          m++;
        }
        for (l=0; l<nset; l++) {
          //printf("setname %s\t vname %s\n",ha_set[l].setname,vname);
          if (strncmp(ha_set[l].setname,vname,m)==0) {
            strcpy(record[i].dimsets[0],vname);
            record[i].dims[0]=ha_set[l].size;
            //printf("size %d\n",record[i].size);
            break;
          }
        }
        readitem = strtok(NULL,")");
        strcpy(vname,readitem);
        m=0;
        while (vname[m] != '\0') {
          m++;
        }
        for (l=0; l<nset; l++) {
          //printf("setname %s\t vname %s\n",ha_set[l].setname,vname);
          if (strncmp(ha_set[l].setname,vname,m)==0) {
            strcpy(record[i].dimsets[1],vname);
            record[i].dims[1]=ha_set[l].size;
            //printf("size %d\n",record[i].size);
            break;
          }
        }
        record[i].begadd=j;
        record[i].size=2;
        record[i].antidims[1]=1;
        record[i].antidims[0]=record[i].dims[1];
        record[i].matsize=record[i].dims[0]*record[i].dims[1];
        j=j+record[i].matsize;
      }
      if (dcount>1) {
        //printf("vname %s\n",vname);
        add=1;
        readitem = strtok(vname,"(");
        strcpy(record[i].cofname,readitem);
        for(k=0; k<dcount; k++) {
          readitem = strtok(NULL,",");
          strcpy(vname,readitem);
          m=0;
          while (vname[m] != '\0') {
            m++;
          }
          for (l=0; l<nset; l++) {
            //printf("setname %s\t vname %s\n",ha_set[l].setname,vname);
            if (strncmp(ha_set[l].setname,vname,m)==0) {
              strcpy(record[i].dimsets[k],vname);
              record[i].dims[k]=ha_set[l].size;
              add=add*ha_set[l].size;
              //printf("size %d\n",record[i].size);
              break;
            }
          }
        }
        readitem = strtok(NULL,")");
        strcpy(vname,readitem);
        m=0;
        while (vname[m] != '\0') {
          m++;
        }
        for (l=0; l<nset; l++) {
          //printf("setname %s\t vname %s\n",ha_set[l].setname,vname);
          if (strncmp(ha_set[l].setname,vname,m)==0) {
            strcpy(record[i].dimsets[k],vname);
            record[i].dims[k]=ha_set[l].size;
            add=add*ha_set[l].size;
            //printf("size %d\n",record[i].size);
            break;
          }
        }
        record[i].begadd=j;
        record[i].size=dcount+1;
        record[i].matsize=add;
        record[i].antidims[dcount]=1;
        for (l=record[i].size-2; l>-1; l--) {
          record[i].antidims[l]=record[i].antidims[l+1]*record[i].dims[l+1];
        }
        j=j+add;
      }
    }
  }*/
  //free(tempcofname);
  //for (i=0;i<ncof;i++) printf("oday %s\n",record[i].cofname);
  return addi;
}

/*int hcge_rsumele(hcge_sumcof *ha_cof,uvadd ncof, ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele,ha_cgesumele *ha_cofele) {
  uvdim l,m,dcount1,dcount,snm;//,l3=0,l4=0,ld2=0,ld3=0,ld4=0,dcountdim[MAXVARDIM]
  uvadd ld1=0,l2=0,l1=0,i,j,sizeele,dcountdim1[MAXVARDIM];
  char vname[NAMESIZE], cofname[NAMESIZE], cofname1[NAMESIZE],interchar[NAMESIZE],sname[NAMESIZE];//vnamecopy[NAMESIZE],, cofname2[NAMESIZE], cofname3[NAMESIZE], cofname4[NAMESIZE]
  //unsigned short int snum;
  //char *readitem=NULL;
  j=0;
  //for (i=0;i<ncof;i++) printf("varname %s\n",ha_cof[i].cofname);
  for (i=0; i<ncof; i++) {
    strcpy(vname,ha_cof[i].sumname);
    //printf("vname %s\n",vname);
    //strcpy(vnamecopy,vname);
    //dcount=ha_cgenchf(vname,')');
    dcount1=ha_cof[i].size;
    if (dcount1==0) {
      strcpy(ha_cofele[j].varname,vname);
      j=j+1;
    } else {
      //dcount=ha_cgenchf(vname,',');
      //dcount=ha_cof[i].size-1;
      //printf("vname %s dcount %d i %d\n",vname,ha_cof[i].size,i);
      if (dcount1==1) {
        //printf("vname %s\n",vname);
        //readitem = strtok(vname,"(");
        //readitem = strtok(NULL,")");
        //strcpy(vname,readitem);
        //m=0;
        //while (vname[m] != '\0') m++;
        for (l=0; l<nset; l++) {
          //printf("setname %s\t vname %s\n",ha_set[l].setname,vname);
          if (strcmp(ha_set[l].setname,ha_cof[i].dimsets[0])==0) {
            //d1=ha_set[l].size;
            //strcpy(cofname,ha_cof[i].cofname);
            strcpy(cofname,ha_cof[i].sumname);
            //readitem = strtok(cofname,"(");
            //strcpy(cofname,readitem);
            strcat(cofname, "(");
            strcpy(cofname1,cofname);
            for (m=0; m<ha_set[l].size; m++) {
              strcpy(cofname,cofname1);
              //sprintf(interchar, "%d", m);
              sprintf(sname, "%d", ha_setele[ha_set[l].begadd+m].setsh);
              //printf("inter %s sname %s\n",interchar,sname);
              strcpy(interchar,sname);
              strcat(cofname, interchar);
              strcat(cofname, ")");
              strcpy(ha_cofele[j+m].varname,cofname);
            }
            j=j+ha_set[l].size;
            break;
          }
        }
      }
      //strcpy(vname,vnamecopy);
      if (dcount1>1) {
        //printf("vname %s\n",vname);
        //readitem = strtok(vname,"(");
        //readitem = strtok(NULL,",");
        //strcpy(vname,readitem);
        //m=0;
        //while (vname[m] != '\0') m++;
        sizeele=1;
        for (dcount=0; dcount<dcount1; dcount++) {
          for (l=0; l<nset; l++) {
            //printf("setname %s\t vname %s\n",ha_cof[i].dimsets[dcount],vname);
            if (strcmp(ha_set[l].setname,ha_cof[i].dimsets[dcount])==0) {
              //d1=ha_set[l].size;
              //dcountdim[dcount]=l;
              sizeele=sizeele*ha_cof[i].dims[dcount];//ha_set[l].size;
              //printf("sizeof %d set %d\n",sizeele,ha_set[l].size);
              break;
            }
          }
        }
        for (dcount=dcount1; dcount>0; dcount--) {
          if(dcount==dcount1) {
            dcountdim1[dcount-1]=1;
          } else {
            dcountdim1[dcount-1]=ha_cof[i].dims[dcount]*dcountdim1[dcount];  //ha_set[dcountdim[dcount]].size
          }
        }
        for (l1=0; l1<sizeele; l1++) {
          strcpy(cofname,ha_cof[i].sumname);
          strcat(cofname, "(");
          //printf("cofele %s\n",cofname);
          //printf("size %d\n",sizeele);
          l2=l1;
          for (dcount=0; dcount<dcount1-1; dcount++) {
            //printf("size %g\n",ha_set[dcountdim[dcount]].begadd+floor(l2/dcountdim1[dcount]));
            snm=(int) l2/dcountdim1[dcount];//ha_set[dcountdim[dcount]].begadd+floor(l2/dcountdim1[dcount]);
            sprintf(sname, "%d",snm);//sprintf(sname, "%d",ha_setele[ld1].setsh);
            //printf("sname %s\n",sname);
            strcpy(interchar,sname);
            //printf("sname %s\n",sname);
            strcat(cofname, interchar);
            strcat(cofname, ",");
            l2=l2-ld1*dcountdim1[dcount];
          }
          snm=(int) l2/dcountdim1[dcount];//ha_set[dcountdim[dcount]].begadd+floor(l2/dcountdim1[dcount]);
          //printf("cofele %s\n",cofname);
          sprintf(sname, "%d",snm);//sprintf(sname, "%d",ha_setele[ld1].setsh);
          strcpy(interchar,sname);
          strcat(cofname, interchar);
          strcat(cofname, ")");
          strcpy(ha_cofele[j+l1].varname,cofname);
        }
        j=j+sizeele;
      }
    }
  }
  return 1;
}*/

int hcge_rcofele(hcge_cof *ha_cof,uvadd ncof, ha_cgeset *ha_set,uvdim nset, ha_cgecofele *ha_cofele) {
  uvdim l,m,dcount;//,dcount1;//,l3=0,l4=0,ld2=0,ld3=0,ld4=0
  uvadd l1=0,i;//,sizeele,dcountdim1[MAXVARDIM],l2=0,dcountdim[MAXVARDIM];//,ld1=0
  //char vname[NAMESIZE], cofname[NAMESIZE], cofname1[NAMESIZE];//,interchar[4];//,sname[NAMESIZE];//vnamecopy[NAMESIZE],, cofname2[NAMESIZE], cofname3[NAMESIZE], cofname4[NAMESIZE]
  //unsigned short int snum;
  //char *readitem=NULL;
  //j=0;
  //for (i=0;i<ncof;i++) printf("varname %s\n",ha_cof[i].cofname);
  for (i=0; i<ncof; i++) {
    for (dcount=0; dcount<ha_cof[i].size; dcount++) {
      for (l=0; l<nset; l++) {
        //printf("setname %s\t vname %s\n",ha_cof[i].dimnames[dcount],vname);
        if (strcmp(ha_set[l].setname,ha_set[ha_cof[i].setid[dcount]].setname)==0) {
          ha_set[ha_cof[i].setid[dcount]].begadd=ha_set[l].begadd;
          //d1=ha_set[l].size;
          //dcountdim[dcount]=l;
          //sizeele=sizeele*ha_set[l].size;
          //printf("sizeof %d set %d\n",sizeele,ha_set[l].size);
          break;
        }
      }
    }
    //strcpy(vname,ha_cof[i].cofname);
    //printf("vname %s\n",vname);
    //strcpy(vnamecopy,vname);
    //dcount=ha_cgenchf(vname,')');
    //dcount1=ha_cof[i].size;
    //if (dcount1==0) {
    //strcpy(ha_cofele[j].cofname,vname);
    //j=j+1;
    //} else {
    //dcount=ha_cgenchf(vname,',');
    //dcount=ha_cof[i].size-1;
    //printf("vname %s dcount %d i %d\n",vname,ha_cof[i].size,i);
    //if (dcount1==1) {
    //printf("vname %s\n",vname);
    //readitem = strtok(vname,"(");
    //readitem = strtok(NULL,")");
    //strcpy(vname,readitem);
    //m=0;
    //while (vname[m] != '\0') m++;
    //for (l=0; l<nset; l++) {
    //printf("setname %s\t vname %s\n",ha_set[l].setname,vname);
    //if (strcmp(ha_set[l].setname,ha_cof[i].dimsets[0])==0) {
    //d1=ha_set[l].size;
    //strcpy(cofname,ha_cof[i].cofname);
    //strcpy(cofname,ha_cof[i].cofname);
    //ha_cof[i].setbegadd1[0]=ha_set[l].begadd;
    //readitem = strtok(cofname,"(");
    //strcpy(cofname,readitem);
    //strcat(cofname, "(");
    //strcpy(cofname1,cofname);
    //for (m=0; m<ha_set[l].size; m++) {
    //strcpy(cofname,cofname1);
    //sprintf(interchar, "%d", m);
    //snum=sprintf(sname, "%d", ha_setele[ha_set[l].begadd+m].setsh);
    //strcpy(interchar,sname);
    //strcat(cofname, interchar);
    //strcat(cofname,",");
    //strcat(cofname, ")");
    //strcpy(ha_cofele[j+m].cofname,cofname);
    //}
    //j=j+ha_set[l].size;
    //break;
    //}
    //}
    //}
    //strcpy(vname,vnamecopy);
    //if (dcount1>1) {
    //printf("vname %s\n",vname);
    //readitem = strtok(vname,"(");
    //readitem = strtok(NULL,",");
    //strcpy(vname,readitem);
    //m=0;
    //while (vname[m] != '\0') m++;
    //sizeele=1;
    //for (dcount=0; dcount<dcount1; dcount++) {
    //for (l=0; l<nset; l++) {
    //printf("setname %s\t vname %s\n",ha_cof[i].dimnames[dcount],vname);
    //if (strcmp(ha_set[l].setname,ha_cof[i].dimsets[dcount])==0) {
    //ha_cof[i].setbegadd1[dcount]=ha_set[l].begadd;
    //d1=ha_set[l].size;
    //dcountdim[dcount]=l;
    //sizeele=sizeele*ha_set[l].size;
    //printf("sizeof %d set %d\n",sizeele,ha_set[l].size);
    //break;
    //}
    //}
    //}
    /*for (dcount=dcount1; dcount>0; dcount--) {
      if(dcount==dcount1) {
        dcountdim1[dcount-1]=1;
      } else {
        dcountdim1[dcount-1]=ha_set[dcountdim[dcount]].size*dcountdim1[dcount];
      }
    }
    for (l1=0; l1<sizeele; l1++) {
      //strcpy(cofname,ha_cof[i].cofname);
      //strcat(cofname, "(");
      //printf("cofele %s\n",cofname);
      //printf("size %d\n",sizeele);
      l2=l1;
      for (dcount=0; dcount<dcount1-1; dcount++) {
        //printf("size %g\n",ha_set[dcountdim[dcount]].begadd+floor(l2/dcountdim1[dcount]));
        //ld1=ha_set[dcountdim[dcount]].begadd+(long int) l2/dcountdim1[dcount];
        //snum=sprintf(sname, "%d",ha_setele[ld1].setsh);
        //strcpy(interchar,sname);
        //strcat(cofname, interchar);
        //strcat(cofname, ",");
        l2=l2-(uvadd) floor(l2/dcountdim1[dcount])*dcountdim1[dcount];
      }
      //ld1=ha_set[dcountdim[dcount]].begadd+(long int) l2/dcountdim1[dcount];
      //snum=sprintf(sname, "%d",ha_setele[ld1].setsh);
      //strcpy(interchar,sname);
      //strcat(cofname, interchar);
      //strcat(cofname, ")");
      //strcpy(ha_cofele[j+l1].cofname,cofname);
    }
    j=j+sizeele;*/
    //}
    //}
  }
  return 1;
}

uvdim ha_cgenset(char *fname) {
  FILE * filehandle;
  char line[TABREADLINE]="\0";
  char *commsyntax="set";
  uvdim j=0;
  filehandle = fopen(fname,"r");
  while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
    j++;
  }
  fclose(filehandle);
  return j;
}

int hcge_rinterset(char *fname, int niodata, hcge_iodata *iodata, ha_cgeset *record,uvdim nset) {
  FILE * filehandle;//, *fileout;
  char line[TABREADLINE]="\0",linecopy[TABREADLINE];//,line1[TABREADLINE];longname[TABREADLINE],
  char header[HEADERSIZE],floginame[NAMESIZE];
  char commsyntax[TABREADLINE],varname[NAMESIZE];
  uvdim j=0,i,intindx[4],sign[2],intvar[2];//,inttype[4]
  int k0,k1=0;
  char *readitem=NULL;//,*p;
  for (j=0; j<nset; j++) {
    if(ha_cgefind(record[j].readele,"intertemporal")>-1) {
      intindx[0]=0;
      intindx[1]=0;
      intindx[2]=0;
      intindx[3]=0;
      //inttype[0]=0;
      //inttype[1]=0;
      //inttype[2]=0;
      //inttype[3]=0;
      sign[0]=0;
      sign[1]=0;
      strcpy(line,record[j].readele);
      //printf("line %s\n",line);
      while (ha_cgefrstr(line," ", ""));
      strcpy(linecopy,line);
      k0=ha_cgenfind(line,",");
      if(k0==1) {
        readitem = strtok(line,",");
        readitem = strtok(NULL,"\0");
        //printf("line %s\n",readitem);
        if(strchr(readitem,'+')!=NULL||strchr(readitem,'-')!=NULL) {
          if (strchr(readitem,'+')!=NULL) {
            sign[0]=1;
          } else {
            sign[0]=2;
          }
          //printf("line %s\n",readitem);
          readitem = strtok(readitem,"-+");
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[0]=atoi(readitem);
          } else {
            //printf("line %s\n",readitem);
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              strcpy(floginame,readitem);
              //printf("flog %s\n",floginame);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            //printf("intindx %d\n",k1);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (ha_cgenchf(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[0]=atoi(readitem);
                    //printf("intindx %d\n",intindx[0]);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
          //printf("read %s\n",line);
          strcpy(line,linecopy);
          //printf("read %s\n",line);
          readitem = strtok(line,",");
          readitem = strtok(NULL,"-+");
          //printf("read %s\n",readitem);
          readitem = strtok(NULL,"\0");
          //printf("read %s\n",readitem);
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[1]=atoi(readitem);
          } else {
            //printf("line %s\n",readitem);
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              strcpy(floginame,readitem);
              //printf("flog %s\n",floginame);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            //printf("intindx %d\n",k1);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (ha_cgenchf(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[1]=atoi(readitem);
                    //printf("intindx %d\n",intindx[0]);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
        } else {
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[0]=atoi(readitem);
          } else {
            //printf("line %s\n",readitem);
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              strcpy(floginame,readitem);
              //printf("flog %s\n",floginame);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            //printf("intindx %d\n",k1);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (ha_cgenchf(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[0]=atoi(readitem);
                    //printf("intindx %d\n",intindx[0]);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
        }
        if(sign[0]>0) if(sign[0]==1) {
            sprintf(line, "%d", intindx[0]+intindx[1]);
            strcpy(record[j].readele,line);
          } else {
            sprintf(line, "%d", intindx[0]-intindx[1]);
            strcpy(record[j].readele,line);
          }
        else {
          sprintf(line, "%d", intindx[0]);
          strcpy(record[j].readele,line);
        }
        record[j].size=1;
      } else {
        //printf("line %s\n",line);
        strcpy(line,linecopy);
        readitem = strtok(line,",");
        readitem = strtok(NULL,",");
        //printf("line %s\n",readitem);
        if(strchr(readitem,'+')!=NULL||strchr(readitem,'-')!=NULL) {
          if (strchr(readitem,'+')!=NULL) {
            sign[0]=1;
          } else {
            sign[0]=2;
          }
          //printf("line %s\n",readitem);
          readitem = strtok(readitem,"-+");
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[0]=atoi(readitem);
          } else {
            //printf("line %s\n",readitem);
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              strcpy(floginame,readitem);
              //printf("flog %s\n",floginame);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            //printf("intindx %d\n",k1);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (ha_cgenchf(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[0]=atoi(readitem);
                    //printf("intindx %d\n",intindx[0]);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
          //printf("read %s\n",line);
          strcpy(line,linecopy);
          //printf("read %s\n",line);
          readitem = strtok(line,",");
          readitem = strtok(NULL,"-+");
          //printf("read %s\n",readitem);
          readitem = strtok(NULL,",");
          //printf("read %s\n",readitem);
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[1]=atoi(readitem);
          } else {
            //printf("line %s\n",readitem);
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              strcpy(floginame,readitem);
              //printf("flog %s\n",floginame);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            //printf("intindx %d\n",k1);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (ha_cgenchf(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[1]=atoi(readitem);
                    //printf("intindx %d\n",intindx[0]);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
        } else {
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[0]=atoi(readitem);
          } else {
            //printf("line %s\n",readitem);
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              strcpy(floginame,readitem);
              //printf("flog %s\n",floginame);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            //printf("intindx %d\n",k1);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (ha_cgenchf(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[0]=atoi(readitem);
                    //printf("intindx %d\n",intindx[0]);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
        }
        if(sign[0]>0) {
          if(sign[0]==1) {
            intvar[0]=intindx[0]+intindx[1];
          } else {
            intvar[0]=intindx[0]+intindx[1];
          }
        } else {
          intvar[0]=intindx[0];
        }
        strcpy(line,linecopy);
        readitem = strtok(line,",");
        readitem = strtok(NULL,",");
        readitem = strtok(NULL,"\0");
        //printf("line1 %s\n",readitem);
        if(strchr(readitem,'+')!=NULL||strchr(readitem,'-')!=NULL) {
          if (strchr(readitem,'+')!=NULL) {
            sign[1]=1;
          } else {
            sign[1]=2;
          }
          //printf("line %s\n",readitem);
          readitem = strtok(readitem,"-+");
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[2]=atoi(readitem);
          } else {
            //printf("line %s\n",readitem);
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              strcpy(floginame,readitem);
              //printf("flog %s\n",floginame);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            //printf("intindx %d\n",k1);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              //printf("datline %s\n",line);
              if (ha_cgenchf(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                //printf("datline %s\n",line);
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[2]=atoi(readitem);
                    //printf("intindx %d\n",intindx[2]);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
          //printf("read %s\n",line);
          strcpy(line,linecopy);
          //printf("read %s\n",line);
          readitem = strtok(line,",");
          readitem = strtok(NULL,",");
          readitem = strtok(NULL,"-+");
          //printf("read %s\n",readitem);
          readitem = strtok(NULL,"\0");
          //printf("read %s\n",readitem);
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[3]=atoi(readitem);
          } else {
            //printf("line %s\n",readitem);
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              strcpy(floginame,readitem);
              //printf("flog %s\n",floginame);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            //printf("intindx %d\n",k1);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (ha_cgenchf(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[3]=atoi(readitem);
                    //printf("intindx %d\n",intindx[0]);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
        } else {
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[2]=atoi(readitem);
          } else {
            //printf("line %s\n",readitem);
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              strcpy(floginame,readitem);
              //printf("flog %s\n",floginame);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            //printf("intindx %d\n",k1);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (ha_cgenchf(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[2]=atoi(readitem);
                    //printf("intindx %d\n",intindx[2]);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
        }
        if(sign[1]>0) if(sign[1]==1) {
            intvar[1]=intindx[2]+intindx[3];
          } else {
            intvar[1]=intindx[2]-intindx[3];
          }
        else {
          intvar[1]=intindx[2];
        }
        record[j].readele[0]='\0';
        for (i=intvar[0]; i<intvar[1]+1; i++) {
          sprintf(line, "%d",i);
          strcat(record[j].readele,line);
          strcat(record[j].readele,",");
        }
        record[j].readele[strlen(record[j].readele)-1]='\0';
        record[j].size=intvar[1]+1-intvar[0];
      }
      //printf("record %s\n",record[j].readele);
    }
  }
  return 0;
}

int ha_cgerset(char *fname, int niodata, hcge_iodata *iodata, ha_cgeset *record,uvdim nset) {
  FILE * filehandle;//, *fileout;
  char line[TABREADLINE]="\0",linecopy[TABREADLINE],line1[TABREADLINE];
  char longname[TABREADLINE],tempvar[256];
  char *commsyntax="set";
  uvdim j=0,vsize,dim1,i;
  int k0,k1,k2;
  char *readitem=NULL;//,*p;

  filehandle = fopen(fname,"r");
  //fileout= fopen("my.txt", "wb");

  while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
    //printf("%s\n",line);
    /*while (ha_cgefrstr(line,"\n", ""));
    while (ha_cgefrstr(line,"\r", ""));
    while (ha_cgefrstr(line,"  ", " "));
    strcat(line, "\n");
    while (ha_cgedrcmt(line,"!"));
    while (ha_cgedrcmt(line,"#"));
    while (ha_cgefrstr(line,"  ", " "));*/
    strcpy(linecopy,line);
    k2=ha_cgefind(line,"intertemporal");
    if(k2>-1) {
      readitem = strtok(line+k2," ");
      readitem = strtok(NULL," ");
      strcpy(record[j].setname,readitem);
      strcpy(line1,"intertemporal,");
      readitem = strtok(NULL,"[");
      readitem = strtok(NULL,"]");
      strcat(line1,readitem);
      readitem = strtok(NULL,";");
      //printf("read %s\n",readitem);
      if (strchr(readitem,'[')!=NULL) {
        strcat(line1,",");
        readitem = strtok(readitem,"[");
        readitem = strtok(NULL,"]");
        strcat(line1,readitem);
      }
      //printf("line1 %s\n",line1);
      strcpy(record[j].readele,line1);
      record[j].intertemp=true;
    } else {
      if (strchr(line,'=')!=NULL) {
        ha_cgefrstr(line,"union", "^");
        readitem = &line[4];
        strcpy(line1,readitem);
        while (ha_cgefrstr(line1," ", ""));
        readitem = strtok(line1,"=");
        //printf("line %s\n",line);
        strcpy(record[j].setname,readitem);
        if (strchr(line,'-')!=NULL) {
          line1[0]='-';
          line1[1]=',';
          line1[2]='\0';
          readitem = strtok(NULL,"-");
          for (i=0; i<nset; i++) {
            if (strcmp(readitem,record[i].setname)==0) {
              dim1=record[i].size;
              strcat(line1,readitem);
              strcat(line1,",");
              break;
            }
          }
          readitem = strtok(NULL,";");
          for (i=0; i<nset; i++) {
            if (strcmp(readitem,record[i].setname)==0) {
              record[j].size=dim1-record[i].size;
              strcat(line1,readitem);
              break;
            }
          }
          strcpy(record[j].readele,line1);
          //printf("%s\n",record[j].readele);
        } else {
          if(strchr(line,'+')!=NULL) {
            line1[0]='+';
            line1[1]=',';
            line1[2]='\0';
            readitem = strtok(NULL,"+");
            for (i=0; i<nset; i++) {
              if (strcmp(readitem,record[i].setname)==0) {
                dim1=record[i].size;
                strcat(line1,readitem);
                strcat(line1,",");
                break;
              }
            }
            readitem = strtok(NULL,";");
            for (i=0; i<nset; i++) {
              if (strcmp(readitem,record[i].setname)==0) {
                record[j].size=dim1+record[i].size;
                strcat(line1,readitem);
                break;
              }
            }
            strcpy(record[j].readele,line1);
          } else {
            if(strchr(line,'^')!=NULL) {
              line1[0]='^';
              line1[1]=',';
              line1[2]='\0';
              readitem = strtok(NULL,"^");
              for (i=0; i<nset; i++) {
                if (strcmp(readitem,record[i].setname)==0) {
                  dim1=record[i].size;
                  strcat(line1,readitem);
                  strcat(line1,",");
                  break;
                }
              }
              readitem = strtok(NULL,";");
              for (i=0; i<nset; i++) {
                if (strcmp(readitem,record[i].setname)==0) {
                  record[j].size=dim1+record[i].size;
                  strcat(line1,readitem);
                  break;
                }
              }
              strcpy(record[j].readele,line1);
            } else {
              line1[0]='=';
              line1[1]='\0';
              readitem = strtok(NULL,";");
              //printf("read 1%s1\n",readitem);
              for (i=0; i<nset; i++) {
                if (strcmp(readitem,record[i].setname)==0) {
                  record[j].size=record[i].size;
                  strcpy(record[j].header,record[i].header);
                  sprintf(tempvar, "%d",i);
                  strcat(line1,tempvar);
                  //strcat(line1,",");
                  //printf("header %s setname %s s %d\n",record[i].header,record[i].setname,record[i].size);
                  break;
                }
              }
              //printf("line1111 %s\n",line1);
              strcpy(record[j].readele,line1);
            }
          }
        }
      } else {
        //while (ha_cgefrstr(line,"  ", " "));
        k0=ha_cgefind(line,"read elements from file");
        //printf("line %s\n",line);
        if(k0>-1) {
          k1=ha_cgefind(line+k0+24," ");
          strncpy(line1,line+k0+24,k1);
          line1[k1]='\0';
          //printf("line1 %s\n",line1);
          for (k0=0; k0<niodata; k0++) if (strcmp(line1,iodata[k0].logname)==0) {
              break;
            }
          record[j].fileid=k0;
        }
        readitem = strtok(line," ");
        readitem = strtok(NULL," ");
        strcpy(record[j].setname,readitem);
        readitem = strtok(NULL,"\"");
        readitem = strtok(NULL,"\"");
        if (readitem!=NULL) {
          strcpy(record[j].header,readitem);
          record[j].readele[0]='\0';
          ha_cgerdvar1(readitem,iodata[record[j].fileid].filname,&vsize,longname,&dim1);
          record[j].size=dim1;
        } else {
          record[j].header[0]='\0';
          readitem = strtok(linecopy,"(");
          readitem = strtok(NULL,")");
          strcpy(record[j].readele,readitem);
          dim1=ha_cgenchf(readitem, ',');
          record[j].size=dim1+1;
        }
      }
    }
    j++;
  }
  //fwrite(&record,sizeof(record),j,fileout);
  //for (i=0;i<j;i++) printf("%s\t%s\t%s\t%d\n",record[i].fileread,record[i].header,record[i].setname,record[i].size);
  //fclose(fileout);
  fclose(filehandle);
  return 0;
}

uvdim ha_setunion(ha_cgesetele *ha_setele, ha_cgeset *ha_set,uvdim nset,uvdim i) {
  uvdim j,l,n,m,dim1=0,dim2=0,j1,sup1,sup2;
  char line[TABREADLINE],*readitem;
  strcpy(line,ha_set[i].readele);
  //ha_set[i].subsetid=1;
  //printf("read %s\n",line);
  readitem = strtok(line,",");
  readitem = strtok(NULL,",");
  for (l=0; l<nset; l++) {
    if (strcmp(readitem,ha_set[l].setname)==0) {
      dim1=ha_set[l].size;
      for (sup1=1; sup1<MAXSUPSET; sup1++)if(ha_set[l].subsetid[sup1]==-1) {
          ha_set[l].subsetid[sup1]=i;
          break;
        }
      //ha_set[i].supersetsize=dim1;
      break;
    }
  }
  m=0;
  for (n=0; n<dim1; n++) {
    strcpy(ha_setele[ha_set[i].begadd+m].setele,ha_setele[ha_set[l].begadd+n].setele);
    ha_setele[ha_set[i].begadd+m].setsh[0]=n;
    ha_setele[ha_set[l].begadd+m].setsh[sup1]=n;
    m++;
  }
  readitem = strtok(NULL,",");
  //printf("read %s\n",readitem);
  for (j=0; j<nset; j++) {
    if (strcmp(readitem,ha_set[j].setname)==0) {
      dim2=ha_set[j].size;
      for (sup2=1; sup2<MAXSUPSET; sup2++)if(ha_set[j].subsetid[sup2]==-1) {
          ha_set[j].subsetid[sup2]=i;
          //printf("read union %s set %s sup2 %d i %d\n",line,ha_set[j].setname,sup2,i);
          break;
        }
      //ha_set[j].subsetid=1;
      break;
    }
  }
  for (n=0; n<dim2; n++) {
    for (j1=0; j1<dim1; j1++) if(strcmp(ha_setele[ha_set[l].begadd+j1].setele,ha_setele[ha_set[j].begadd+n].setele)==0) break;
    if(j1==dim1) {
      strcpy(ha_setele[ha_set[i].begadd+m].setele,ha_setele[ha_set[j].begadd+n].setele);
      ha_setele[ha_set[i].begadd+m].setsh[0]=m;
      ha_setele[ha_set[j].begadd+n].setsh[sup2]=m;
      //ha_setele[ha_set[j].begadd+n].setsh=m;
      //printf("set %s begadd %d sup2 %d i %d m %d n %d in sup %d\n",ha_set[j].setname,ha_set[j].begadd+n,sup2,i,m,n,ha_setele[ha_set[j].begadd+n].setsh[sup2]);
      m++;
    } else {
      ha_setele[ha_set[j].begadd+n].setsh[sup2]=j1;
    }
  }
  ha_set[i].size=m;
  return m;
}
uvdim ha_setplus(ha_cgesetele *ha_setele, ha_cgeset *ha_set,uvdim nset,uvdim i) {
  uvdim j,l,n,m,dim1=0,dim2=0,j1,sup1,sup2;
  char line[TABREADLINE],*readitem;
  strcpy(line,ha_set[i].readele);
  //ha_set[i].subsetid=1;
  //printf("read %s\n",line);
  readitem = strtok(line,",");
  readitem = strtok(NULL,",");
  for (l=0; l<nset; l++) {
    if (strcmp(readitem,ha_set[l].setname)==0) {
      dim1=ha_set[l].size;
      for (sup1=1; sup1<MAXSUPSET; sup1++)if(ha_set[l].subsetid[sup1]==-1) {
          ha_set[l].subsetid[sup1]=i;
          break;
        }
      //ha_set[i].supersetsize=dim1;
      break;
    }
  }
  m=0;
  for (n=0; n<dim1; n++) {
    strcpy(ha_setele[ha_set[i].begadd+m].setele,ha_setele[ha_set[l].begadd+n].setele);
    ha_setele[ha_set[i].begadd+m].setsh[0]=n;
    ha_setele[ha_set[l].begadd+m].setsh[sup1]=n;
    m++;
  }
  readitem = strtok(NULL,",");
  //printf("read %s\n",readitem);
  for (j=0; j<nset; j++) {
    if (strcmp(readitem,ha_set[j].setname)==0) {
      dim2=ha_set[j].size;
      for (sup2=1; sup2<MAXSUPSET; sup2++)if(ha_set[j].subsetid[sup2]==-1) {
          ha_set[j].subsetid[sup2]=i;
          break;
        }
      //ha_set[j].subsetid=1;
      break;
    }
  }
  for (n=0; n<dim2; n++) {
    strcpy(ha_setele[ha_set[i].begadd+m].setele,ha_setele[ha_set[j].begadd+n].setele);
    ha_setele[ha_set[i].begadd+m].setsh[0]=m;
    ha_setele[ha_set[j].begadd+n].setsh[sup2]=m;
    m++;
  }
  return m;
}
uvdim ha_setminus(ha_cgesetele *ha_setele, ha_cgeset *ha_set,uvdim nset,uvdim i) {
  uvdim j,l,n,m,v,indi,dim1=0,dim2=0,sup1;
  char line[TABREADLINE],*readitem;
  strcpy(line,ha_set[i].readele);
  //ha_set[i].subsetid=1;
  readitem = strtok(line,",");
  readitem = strtok(NULL,",");
  for (l=0; l<nset; l++) {
    if (strcmp(readitem,ha_set[l].setname)==0) {
      dim1=ha_set[l].size;
      //printf("name %s\n",readitem);
      for (sup1=1; sup1<MAXSUPSET; sup1++)if(ha_set[i].subsetid[sup1]==-1) {
          ha_set[i].subsetid[sup1]=l;
          //printf("name %s\n",readitem);
          break;
        }
      //ha_set[i].supersetsize=dim1;
      break;
    }
  }
  readitem = strtok(NULL,",");
  for (j=0; j<nset; j++) {
    if (strcmp(readitem,ha_set[j].setname)==0) {
      dim2=ha_set[j].size;
      //ha_set[i].size=dim1-dim2;
      break;
    }
  }
  m=0;
  for (n=0; n<dim1; n++) {
    indi=0;
    for (v=0; v<dim2; v++) {
      if (strcmp(ha_setele[ha_set[l].begadd+n].setele,ha_setele[ha_set[j].begadd+v].setele)==0) {
        indi=1;
        break;
      }
    }
    if(indi==0) {
      strcpy(ha_setele[ha_set[i].begadd+m].setele,ha_setele[ha_set[l].begadd+n].setele);
      ha_setele[ha_set[i].begadd+m].setsh[sup1]=n;
      ha_setele[ha_set[i].begadd+m].setsh[0]=m;
      m++;
    }
  }
  //ha_set[i].size=m;
  return m;
}
uvadd ha_cgersubset(char *fname, ha_cgesetele *ha_setele, ha_cgeset *ha_set,uvdim nset) {
  FILE * filehandle;//, *fileout;
  char line[TABREADLINE]="\0";
  char set[NAMESIZE],subset[NAMESIZE];
  char *commsyntax="subset";
  uvdim i,setd,subsetd,sup1;//,nlength;
  uvadd jj,jjj,j=0,succ=0,ssize=0;
  char *readitem=NULL;

  filehandle = fopen(fname,"r");
  //fileout= fopen("my.txt", "wb");

  while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
    //printf("%s\n",line);
    /*while (ha_cgefrstr(line,"\n", ""));
    while (ha_cgefrstr(line,"\r", ""));
    while (ha_cgefrstr(line,"  ", " "));
    strcat(line, "\n");
    while (ha_cgedrcmt(line,"!"));
    while (ha_cgedrcmt(line,"#"));
    while (ha_cgefrstr(line,"  ", " "));*/
    while (ha_cgefrstr(line," ;", ";"));
    readitem = strtok(line," ");
    readitem = strtok(NULL," ");
    strcpy(subset,readitem);
    //subsetd=0;
    //while (subset[subsetd] != '\0') {
    //subsetd++;
    //}
    readitem = strtok(NULL," ");
    readitem = strtok(NULL," ");
    readitem = strtok(NULL," ");
    readitem = strtok(NULL,";");
    strcpy(set,readitem);
    //setd=0;
    //while (set[setd] != '\0') {
    //setd++;
    //}
    //printf("line %s subset %s set %s\n",line,subset,set);
    succ=0;
    for (i=0; i<nset; i++) {
      if (strcmp(ha_set[i].setname,subset)==0) { //,subsetd)==0) {
        for (j=0; j<nset; j++) {
          if (strcmp(ha_set[j].setname,set)==0) { //,setd)==0) {
            ssize=ha_set[j].size;
            for (sup1=1; sup1<MAXSUPSET; sup1++)if(ha_set[i].subsetid[sup1]==-1) {
                ha_set[i].subsetid[sup1]=j;
                break;
              }
            if(sup1==MAXSUPSET)printf("Error!!! Number of superset exceeded! Increase MAXSUPSET!\n");
            //ha_set[i].subsetid=1;
            //ha_set[i].supersetsize=ha_set[j].size;
            for (jj=ha_set[j].begadd; jj<ha_set[j].begadd+ha_set[j].size; jj++) {
              for (jjj=ha_set[i].begadd; jjj<ha_set[i].begadd+ha_set[i].size; jjj++) {
                //nlength=0;
                //while (ha_setele[jjj].setele[nlength] != '\0') nlength++;
                if (strcmp(ha_setele[jj].setele,ha_setele[jjj].setele)==0) {
                  //printf("set %s\n subset %s\n",ha_setele[jj].setele,ha_setele[jjj].setele);
                  //printf("jjj %d jj %d\n",jjj,jj);
                  //l=sprintf(sname, "%d", jj);
                  //strcpy(ha_setele[jjj].setsh,sname);
                  //ha_setele[jjj].setsh[l]='\0';
                  ha_setele[jjj].setsh[sup1]=ha_setele[jj].setsh[0];
                  succ++;
                  break;
                }
              }
            }
            break;
          }
        }
        if(j==nset)printf("Error!!! No such set %s\n",set);
        break;
      }
      j++;
    }
    if(i==nset)printf("Error!!! No such set %s\n",subset);
    if(succ-ha_set[i].size!=0)printf("Error!!! No set or some elements of set %s are not elements of set %s\n",subset,set);
    // fwrite(&record,sizeof(record),j,fileout);
    // for (i=0;i<j;i++) printf("%s\t%s\t%s\t%d\n",record[i].fileread,record[i].header,record[i].setname,record[i].size);
    // fclose(fileout);
  }
  fclose(filehandle);
  return j;
}

uvdim ha_cgesubsetchck(ha_cgesetele *ha_setele, ha_cgeset *ha_set,uvdim nset,uvadd* contin) {
  uvdim i,setd,subsetd,sup1,sup2,sup3,b;//,nlength;
  uvadd jj,jjj,j=0;
  *contin=0;
  for (i=0; i<nset; i++) {
    b=1;
    while(b==1) {
      for (sup1=1; sup1<MAXSUPSET; sup1++) {
        b=0;
        if(ha_set[i].subsetid[sup1]==-1) {
          b=2;
          break;
        }
        //if(i==nset-1)printf("beg b %d id %d jjj %ld\n",b,ha_set[i].subsetid[sup1],jjj);
        for (sup2=1; sup2<MAXSUPSET; sup2++) {
          b=0;
          if(ha_set[ha_set[i].subsetid[sup1]].subsetid[sup2]==-1) {
            b=2;
            break;
          }
        //if(i==nset-1)printf("beg b1 %d\n",b);
          for (sup3=1; sup3<MAXSUPSET; sup3++){
          //if(i==nset-1)printf("b %d sup1 %d sup2 %d sup3 %d subid1 %d subid2 %d\n",b,sup1,sup2,sup3,ha_set[ha_set[i].subsetid[sup1]].subsetid[sup2],ha_set[i].subsetid[sup3]);
            if(ha_set[ha_set[i].subsetid[sup1]].subsetid[sup2]==ha_set[i].subsetid[sup3]) {
              b=2;
              break;
            }
          }
          //if(sup3==MAXSUPSET)b=0;
          if(b==0)break;
        }
        if(b==0)break;
      }
      if(b==0) {
        *contin=1;
        j=ha_set[ha_set[i].subsetid[sup1]].subsetid[sup2];
        for (sup1=1; sup1<MAXSUPSET; sup1++){
          //printf("i %d idx %d subsetid %d\n",i,sup1,ha_set[i].subsetid[sup1]);
          if(ha_set[i].subsetid[sup1]==-1)break;
        }
        if(sup1==MAXSUPSET){
          printf("Error!!!! Superset size exceeded!!\n");
          *contin=0;
          return 0;
        }
        ha_set[i].subsetid[sup1]=j;
        for (jj=ha_set[j].begadd; jj<ha_set[j].begadd+ha_set[j].size; jj++) {
          for (jjj=ha_set[i].begadd; jjj<ha_set[i].begadd+ha_set[i].size; jjj++) {
            if (strcmp(ha_setele[jj].setele,ha_setele[jjj].setele)==0) {
              ha_setele[jjj].setsh[sup1]=ha_setele[jj].setsh[0];
              break;
            }
          }
        }
        //printf("sub set %s sup set %s\n",ha_set[i].setname,ha_set[j].setname);
        b=1;
      }
    }
  }
  for (i=nset-1; i>-1; i--) {
    b=1;
    while(b==1) {
      for (sup1=1; sup1<MAXSUPSET; sup1++) {
        b=0;
        if(ha_set[i].subsetid[sup1]==-1) {
          b=2;
          break;
        }
        //if(i==nset-1)printf("beg b %d id %d jjj %ld\n",b,ha_set[i].subsetid[sup1],jjj);
        for (sup2=1; sup2<MAXSUPSET; sup2++) {
          b=0;
          if(ha_set[ha_set[i].subsetid[sup1]].subsetid[sup2]==-1) {
            b=2;
            break;
          }
        //if(i==nset-1)printf("beg b1 %d\n",b);
          for (sup3=1; sup3<MAXSUPSET; sup3++){
          //if(i==nset-1)printf("b %d sup1 %d sup2 %d sup3 %d subid1 %d subid2 %d\n",b,sup1,sup2,sup3,ha_set[ha_set[i].subsetid[sup1]].subsetid[sup2],ha_set[i].subsetid[sup3]);
            if(ha_set[ha_set[i].subsetid[sup1]].subsetid[sup2]==ha_set[i].subsetid[sup3]) {
              b=2;
              break;
            }
          }
          //if(sup3==MAXSUPSET)b=0;
          if(b==0)break;
        }
        if(b==0)break;
      }
      if(b==0) {
        *contin=1;
        j=ha_set[ha_set[i].subsetid[sup1]].subsetid[sup2];
        for (sup1=1; sup1<MAXSUPSET; sup1++){
          //printf("i %d idx %d subsetid %d\n",i,sup1,ha_set[i].subsetid[sup1]);
          if(ha_set[i].subsetid[sup1]==-1)break;
        }
        if(sup1==MAXSUPSET){
          printf("Error!!!! Superset size exceeded!!\n");
          *contin=0;
          return 0;
        }
        ha_set[i].subsetid[sup1]=j;
        for (jj=ha_set[j].begadd; jj<ha_set[j].begadd+ha_set[j].size; jj++) {
          for (jjj=ha_set[i].begadd; jjj<ha_set[i].begadd+ha_set[i].size; jjj++) {
            if (strcmp(ha_setele[jj].setele,ha_setele[jjj].setele)==0) {
              ha_setele[jjj].setsh[sup1]=ha_setele[jj].setsh[0];
              break;
            }
          }
        }
        //printf("sub set %s sup set %s\n",ha_set[i].setname,ha_set[j].setname);
        b=1;
      }
    }
  }
  return 1;
}
/*
char *ha_cgerdatf(char *commsyntax, FILE *filehandle, char *readline)
{
  int check1=0,i,count1=0;
  while (commsyntax[count1] != '\0') {
    count1++;
  }
  char uppercomsyn[TABLINESIZE],*n;
  char line[DATREADLINE],linecomm[DATREADLINE],upperlinecomm[DATREADLINE],*finditem=";";//,linecopy[TABLINESIZE+2]
  for (i=0; i<count1; i++) {
    uppercomsyn[ i ] = toupper( (int) commsyntax[ i ] );
    //printf("comsyn %d i %d\n",count1,i);
  }
  //printf("before %s\n",readline);
  strcpy(readline,"\0");

  while (fgets(line,DATREADLINE,filehandle)) {
    //printf("before %s\n",line);
    while (ha_cgefrstr(line,"  ", " "));
    if (line[0]==' ') {
      memmove(line,&line[0]+1,strlen(line)-1);
    }
    //strcpy(linecopy,line);
    strncpy(linecomm,line,count1);
    for ( i = 0; i<count1; i++) {
      upperlinecomm[ i ] = toupper( (int) linecomm[ i ] );
    }
    if (strncmp(upperlinecomm,uppercomsyn,count1) == 0&&check1==0) {
      check1=1;
      n=strstr(line,finditem);//ha_cgefendofc
      if (n==NULL) {
        strcpy(readline,line);
      } else {
        strcpy(readline,line);
        return readline;
      }
      while (fgets(line,DATREADLINE,filehandle)) {
        //strcpy(linecopy,line);
        n=strstr(line,finditem);//ha_cgefendofc
        if (n==NULL) {
          strcat(readline, line);
        } else {
          strcat(readline, line);
          return readline;
        }
      }

    }
  }
  return NULL;
}*/

char *ha_cgercls(char *commsyntax, FILE *filehandle, char *readline) {
  int check1=0,i,count1=0;
  while (commsyntax[count1] != '\0') {
    count1++;
  }
  char uppercomsyn[TABLINESIZE],*n,*p;
  char line[TABLINESIZE],linecomm[TABLINESIZE],upperlinecomm[TABLINESIZE],*finditem=";";//,linecopy[TABLINESIZE+2]
  for (i=0; i<count1; i++) {
    uppercomsyn[ i ] = toupper( (int) commsyntax[ i ] );
    //printf("comsyn %d i %d\n",count1,i);
  }
  //printf("before %s\n",readline);
  strcpy(readline,"\0");

  while (fgets(line,TABLINESIZE,filehandle)) {
    //printf("before %s\n",line);
    while (ha_cgefrstr(line,"  ", " "));
    if (line[0]==' ') {
      memmove(line,&line[0]+1,strlen(line)-1);
    }
    //strcpy(linecopy,line);
    strncpy(linecomm,line,count1);
    for ( i = 0; i<count1; i++) {
      upperlinecomm[ i ] = toupper( (int) linecomm[ i ] );
    }
    if (strncmp(upperlinecomm,uppercomsyn,count1) == 0&&check1==0) {
      //printf("before %s\n",line);
      p=strpbrk(line,"!");
      while(p!=NULL) {
        line[p-line]='\0';
        p=strpbrk(line,"!");
      }
      //printf("before1 %s\n",line);
      //p=strpbrk(line,"\n");
      //if (p!=NULL) line[p-line]='\0';
      //ha_cgefrstr(line,"\n", " ");
      //printf("before1 %s\n",line);
      check1=1;
      n=strstr(line,finditem);//ha_cgefendofc
      if (n==NULL) {
        strcpy(readline,line);
      } else {
        strcpy(readline,line);
        //printf("%s\n",readline);
        return readline;
      }
      while (fgets(line,TABLINESIZE,filehandle)) {
        //printf("before2 %s\n",line);
        //strcpy(linecopy,line);
        p=strpbrk(line,"!");
        while(p!=NULL) {
          line[p-line]='\0';
          p=strpbrk(line,"!");
        }
        //printf("before3 %s\n",line);
        //p=strpbrk(line,"\n");
        //if (p!=NULL) line[p-line]='\0';
        //ha_cgefrstr(line,"\n", " ");
        //printf("before3 %s\n",line);
        n=strstr(line,finditem);//ha_cgefendofc
        if (n==NULL) {
          strcat(readline, line);
          //printf("before4 %s\n",readline);
        } else {
          strcat(readline, line);
          //printf("%s\n",readline);
          return readline;
        }
      }

    }
  }
  return NULL;
}

char *ha_cgertabl(char *commsyntax, FILE *filehandle, char *readline,uvadd rlinesize) {
  int check1=0,i,count1=0;
  while (commsyntax[count1] != '\0') {
    count1++;
  }
  char uppercomsyn[TABLINESIZE],*n;
  char line[TABLINESIZE],linecomm[TABLINESIZE],upperlinecomm[TABLINESIZE],*finditem=";";//,linecopy[TABLINESIZE+2]
  for (i=0; i<count1; i++) {
    uppercomsyn[ i ] = toupper( (int) commsyntax[ i ] );
    //printf("comsyn %d i %d\n",count1,i);
  }
  //printf("before %s\n",readline);
  strcpy(readline,"\0");

  while (fgets(line,TABLINESIZE,filehandle)) {
    //printf("before %s\n",line);
    while (ha_cgefrstr(line,"  ", " "));
    if (line[0]==' ') {
      memmove(line,&line[0]+1,strlen(line)-1);
    }
    //strcpy(linecopy,line);
    strncpy(linecomm,line,count1);
    for ( i = 0; i<count1; i++) {
      upperlinecomm[ i ] = toupper( (int) linecomm[ i ] );
    }
    if (strncmp(upperlinecomm,uppercomsyn,count1) == 0&&check1==0) {
      check1=1;
      n=strstr(line,finditem);//ha_cgefendofc
      if (n==NULL) {
        strcpy(readline,line);
      } else {
        strcpy(readline,line);
          if(strlen(readline)>=((size_t)(9*rlinesize/10)))printf("Warning!!!! Overflow may happen in function ha_cgertabl. Increase TABREADLINE or DATREADLINE may help!%ld %ld\n",strlen(readline),rlinesize);
        return readline;
      }
      while (fgets(line,TABLINESIZE,filehandle)) {
        //strcpy(linecopy,line);
        n=strstr(line,finditem);//ha_cgefendofc
        if (n==NULL) {
          strcat(readline, line);
        } else {
          strcat(readline, line);
          if(strlen(readline)>=((size_t)(9*rlinesize/10)))printf("Warning!!!! Overflow may happen in function ha_cgertabl. Increase TABREADLINE or DATREADLINE may help! %ld %ld\n",strlen(readline),rlinesize);
          return readline;
        }
      }

    }
  }
  return NULL;
}
char *ha_cgertabl1(char *commsyntax, FILE *filehandle, char *readline, ha_cgevar *record, hcge_cof *ha_cof,uvadd ncof,ha_cgetype *zerodivide,uvadd rlinesize) {
  int check1=0,count1=0;
  char *zerosyntax="zerodivide default",*p,*zerosyntax1="zerodivide (",*zerosyntax2="zerodivide(";
  while (commsyntax[count1] != '\0') {
    count1++;
  }
  char *n,line[TABLINESIZE],*finditem=";";//,linecopy[TABLINESIZE+2]
  //for ( i = 0; i<count1; i++) {
  //uppercomsyn[ i ] = toupper( (int) commsyntax[ i ] );
  //}
  //printf("before %s\n",readline);
  strcpy(readline,"\0");
  while (fgets(line,TABLINESIZE,filehandle)) {
    while (ha_cgefrstr(line,"  ", " "));
    if (line[0]==' ') {
      memmove(line,&line[0]+1,strlen(line)-1);
    }
    //printf("line %s\n",line);
    //strcpy(linecopy,line);
    //strncpy(linecomm,line,count1);
    //for ( i = 0; i<count1; i++) {
    //upperlinecomm[ i ] = toupper( (int) linecomm[ i ] );
    //}
    if (strncmp(line,commsyntax,count1) == 0&&check1==0) {
      check1=1;
      n=strstr(line,finditem);//ha_cgefendofc
      if (n==NULL) {
        strcpy(readline,line);
      } else {
        strcpy(readline,line);
          if(strlen(readline)>=((size_t)(9*rlinesize/10)))printf("Warning!!!! Overflow may happen in function ha_cgertabl1. Increase TABREADLINE or DATREADLINE may help! %ld %ld\n",strlen(readline),rlinesize);
        return readline;
      }
      while (fgets(line,TABLINESIZE,filehandle)) {
        //strcpy(linecopy,line);
        n=strstr(line,finditem);//ha_cgefendofc
        if (n==NULL) {
          strcat(readline, line);
        } else {
          strcat(readline, line);
          if(strlen(readline)>=((size_t)(9*rlinesize/10)))printf("Warning!!!! Overflow may happen in function ha_cgertabl1. Increase TABREADLINE or DATREADLINE may help! %ld %ld\n",strlen(readline),rlinesize);
          return readline;
        }
      }

    }
    if (strncmp(line,commsyntax,count1) != 0) {
      //printf("line %s\n",line);
      //strncpy(linecomm,line,18);
      //for ( i = 0; linecomm[ i ]; i++) {
      //upperlinecomm[ i ] = toupper( (int) linecomm[ i ] );
      //}
      if (strncmp(line,zerosyntax,18) == 0) {
        p=strtok(&line[18],";");
        //printf("p %s\n",p);
        *zerodivide=hnew_simplrpl(p,record,ha_cof,ncof);//atof(p);
        //printf("p %s zero %f\n",p,*zerodivide);
      }
      if (strncmp(line,zerosyntax1,12) == 0||strncmp(line,zerosyntax2,11) == 0) {
        //printf("line %s\n",line);
        p=strtok(line,")");
        p=strtok(NULL,"default");
        p=strtok(NULL," ");
        p=strtok(NULL,";");
        //printf("p %s\n",p);
        *zerodivide=hnew_simplrpl(p,record,ha_cof,ncof);//atof(p);
        //printf("p %s zero %f\n",p,*zerodivide);
      }
      //return NULL;
    }
  }
  return NULL;
}

char *ha_cgefrstr1(char *line, char *finditem, char *replitem) {
  char buffer[DATREADLINE];
  int count2 = 0,index=0;
  while (finditem[count2] != '\0') {
    count2++;
  }
  index=ha_cgefind(line, finditem);
  if (index==-1) {
    return NULL;
  }
  strncpy(buffer, line, index);
  buffer[index] = '\0';
  strcat(buffer, replitem);
  strcat(buffer, line+index+count2);
  strcpy(line,buffer);
  return line;
}

char *ha_cgefrstrvbz1(char *line, char *finditem, char *replitem,uvdim nbuffer) {
  char *buffer= (char *) calloc (nbuffer,sizeof(char));
  //char buffer[TABREADLINE];
  int count2 = 0,index=0;
  while (finditem[count2] != '\0') {
    count2++;
  }
  index=ha_cgefind(line, finditem);
  if (index==-1) {
    free(buffer);
    return NULL;
  }
  strncpy(buffer, line, index);
  buffer[index] = '\0';
  strcat(buffer, replitem);
  strcat(buffer, line+index+count2);
  strcpy(line,buffer);
  free(buffer);
  return line;
}

char *ha_cgefrstr(char *line, char *finditem, char *replitem) {
  char buffer[DATREADLINE];
  char *p;
  unsigned short int count2 = 0,index;
  while (finditem[count2] != '\0') {
    count2++;
  }

  if (!(p = strstr(line, finditem))) {
    return NULL;
  }
  index=p-line;
  strncpy(buffer, line, index);
  buffer[index] = '\0';
  strcat(buffer, replitem);
  strcat(buffer, (line+index+count2));
  strcpy(line,buffer);
  return line;
}

char *ha_cgefrstrvbz(char *line, char *finditem, char *replitem,uvdim nbuffer) {
  char *buffer= (char *) calloc (nbuffer,sizeof(char));
  //char buffer[TABREADLINE];
  char *p;
  long int count2 = 0,index;
  while (finditem[count2] != '\0') {
    count2++;
  }

  if (!(p = strstr(line, finditem))) {
    free(buffer);
    return NULL;
  }
  index=p-line;
  strncpy(buffer, line, index);
  buffer[index] = '\0';
  strcat(buffer, replitem);
  strcat(buffer, (line+index+count2));
  strcpy(line,buffer);
  free(buffer);
  return line;
}

char *ha_cgefrchr(char *line, int finditem, int replitem) {
  char *p;
  if (!(p = strchr(line, finditem))) {
    return NULL;
  }
  *p=replitem;
  return line;
}

int ha_cgefind(char *line, char *finditem) {
  int i;
  char *t1,*t2;
  for (i=0; line[i]; i++) {
    t1 = &line[i];
    t2=&finditem[0];
    while (*t2 && tolower((int)*t2)==tolower((int)*t1)) {
      t1++;
      t2++;
    }
    if (!*t2) {
      return i;
    }
  }
  return -1;
}//ref: http://www.java2s.com/Code/C/String/Findsubstringourownfunction.htm

char *ha_cgedrcmt(char *line, char *token) {
  char buffer[TABREADLINE];
  char *p,*readitem=NULL;
  int count1 = 0;
  while (line[count1] != '\0') {
    count1++;
  }

  if (!(p = strstr(line, token))) {
    return NULL;
  }
  strncpy(buffer, line, p-line);
  buffer[p-line] = '\0';
  strcat(buffer," ");
  //if (p-line==0) printf("line %s\n",line);
  readitem = strtok(line,token);
  //if (p-line==0) printf("line %s\n",readitem);
  if ((p-line)!=0) {
    readitem = strtok(NULL,token);
  }
  //if (token='%') printf("line %s\n",readitem);
  readitem = strtok(NULL,"\0");
  if(readitem!=NULL) {
    strcat(buffer,readitem);
  }
  strcpy(line,buffer);
  return line;
}

/*int ha_spadd_size(int *irnA, int *jcnA, int *irnB, int *jcnB,int nrow,int *irnC){
  int i,i1,i2,j,nrowplus=nrow+1,rnzA,rnzB;
  for(i=0;i<nrow;++i){
    rnzA=irnA[i+1]-irnA[i];
    rnzB=irnB[i+1]-irnB[i];
    if(rnzA==0||rnzB==0){
      irnC[i+1]=irnC[i]+rnzA+rnzB;
      continue;
    }
    if(jcnA[irnA[i+1]-1]<jcnB[irnB[i+1]-1]){
      i1=0;
      i2=0;
        for(j=0; j<rnzA; ++j) {
          while (jcnA[j]>jcnB[i2])++i2;
          if(jcnA[j]<jcnB[i2])++i1;
        }
        irnC[i+1]=irnC[i]+rnzA+i1;
    }else{
      i1=0;
      i2=0;
        for(j=0; j<rnzB; ++j) {
          while (jcnB[j]>jcnA[i2])++i2;
          if(jcnB[j]<jcnA[i2])++i1;
        }
        irnC[i+1]=irnC[i]+rnzB+i1;
    }
  }
  return 1;
}

int ha_spadd_oper(int *irnA, int *jcnA, ha_cgetype *valsA, int *irnB, int *jcnB, ha_cgetype *valsB,int nrow,int *irnC, int *jcnC,ha_cgetype *valsC){
  int i,i1,i2,j,nrowplus=nrow+1,rnzA,rnzB,rnzC,itmp;
  for(i=0;i<nrow;++i){
    rnzA=irnA[i+1]-irnA[i];
    rnzB=irnB[i+1]-irnB[i];
    rnzC=irnC[i+1]-irnC[i];
    if(rnzA==0){
      memcpy(valsC+irnC[i],valsB+irnB[i],rnzB*sizeof(ha_cgetype));
      continue;
    }
    if(rnzB==0){
      memcpy(valsC+irnC[i],valsA+irnA[i],rnzA*sizeof(ha_cgetype));
      continue;
    }
      i1=rnzA-1;
      i2=rnzB-1;
      j=rnzC-1;
      while(j>-1) {
        if(jcnA[i1]==jcnB[i2]) {
          valsC[j]=valsA[i1]+valsB[i2];
          jcnC[j]=jcnA[i1];
          --i1;
          --i2;
          --j;
          continue;
        }
        itmp=i1;
        while(irnA[i1]>irnB[i2]) {
          //vecbivi[i]=vecbivi[j1];
          //biviindx1[i]=biviindx1[j1];
          //i1;
          --i1;
          --j;
        }
        if(i1<itmp){
          memcpy(valsC+irnC[i],valsA+irnA[itmp],(itmp-i1)*sizeof(ha_cgetype));
          memcpy(jcnC+irnC[i],jcnA+irnA[itmp],(itmp-i1)*sizeof(int));
        }
        itmp=i2;
        while(irnB[i2]<irnA[i1]) {
          //vecbivi[i]=vecbivi0[j];
          //biviindx1[i]=biviindx0[j];
          --i2;
          --j;
        }
        if(i1<itmp){
          memcpy(valsC+irnC[i],valsB+irnB[itmp],(itmp-i2)*sizeof(ha_cgetype));
          memcpy(jcnC+irnC[i],jcnB+irnB[itmp],(itmp-i2)*sizeof(int));
        }
      }
  }
  return 1;
}*/
