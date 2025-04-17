//#include <ha_cgefparse.h>
//#include <ha_cgetab.h>
#include <ha_cgeglobal.h>

ha_cgetype hnew_simplrpl(char *var2, ha_cgevar *record, hcge_cof *ha_cof,uvadd ncof) {
  uvadd index;
  while (ha_cgefrstr(var2," ", ""));
  ha_cgetype eval=0;
  if (var2[0]>='0'&&var2[0]<='9') {
    eval=atof(var2);
    return eval;
  }
  index=ncof-1;
  do {
    if (strcmp(ha_cof[index].cofname,var2)==0) {
      eval=record[ha_cof[index].begadd].varval;
      break;
    }
  } while (index--);
  return eval;
}


int hnew_varrepl(char *var2, ha_cgeset *ha_set,hcge_cof *ha_cof,uvadd ncof, hcge_cof *ha_var,uvadd nvar,uvadd ncofele,hcge_sumcof *sum_cof,int totalsum,hcge_calvars *ha_calvar,int ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim,int varindex) {
  uvadd index;
  char *p=NULL;//,copyvar[TABREADLINE];//,*p1=NULL,*p2=NULL,*p3=NULL,*p4=NULL;
  uvdim l1,l,sup;//,svar2;//=0,i2=0,i3=0,i4=0,svar1,svar2,checkvar20=0,checkvar10=0,checkvar11=0,checkvar12=0,checkvar16=0,checkvar14=0,l;
  int leadlag;
  bool IsChange=false;
  //printf("var2 %s\n",var2);
  p= strtok(var2,"{");
  if (p==NULL) {
    p=&var2[0];
  }
  if(var2[0]=='p'&&var2[1]=='_') {
    IsChange=true;
    p++;
    p++;
  }
  //printf("p1 %s\n",p);
  index=ncof-1;
  do {
    if (strcmp(ha_cof[index].cofname,p)==0) {
      if(!ha_cof[index].suplval)printf("Warning!!!! coefficient %s has not been supplied with values!\n",ha_cof[index].cofname);
      if(varindex==2) {
        ha_calvar[ha_calvarsize].Var2BegAdd=ha_cof[index].begadd;
        ha_calvar[ha_calvarsize].Var2Type=0;
      } else {
        ha_calvar[ha_calvarsize].Var1BegAdd=ha_cof[index].begadd;
        ha_calvar[ha_calvarsize].Var1Type=0;
      }
      //printf("p %s cof %s index %d add %d size %d varindex %d\n",p,ha_cof[index].cofname,ha_cof[index].begadd,ha_calvar[ha_calvarsize].Var1BegAdd,ha_calvarsize,varindex);
      switch(ha_cof[index].size) {
      case 0:
        for (l=0; l<fdim; l++) {
          if(varindex==2) {
            ha_calvar[ha_calvarsize].Var2ADims[l]=0;
            //ha_calvar[ha_calvarsize].Var2leadlag[l]=0;
            //ha_calvar[ha_calvarsize].Var1Type=1;
            //ha_calvar[ha_calvarsize].Var1BegAdd=ha_cof[index].begadd;
            //ha_calvar[ha_calvarsize].Var1Size=0;
          } else {
            ha_calvar[ha_calvarsize].Var1ADims[l]=0;
            //ha_calvar[ha_calvarsize].Var1leadlag[l]=0;
            //ha_calvar[ha_calvarsize].Var2Type=1;
            //ha_calvar[ha_calvarsize].Var2BegAdd=ha_cof[index].begadd;
            //ha_calvar[ha_calvarsize].Var2Size=0;
          }
        }
        break;
      case 1:
        //if(varindex==2) ha_calvar[ha_calvarsize].Var2BegAdd=ha_cof[index].begadd;
        //else ha_calvar[ha_calvarsize].Var1BegAdd=ha_cof[index].begadd;
        p=strtok(NULL,"}");
        leadlag=0;
        hnew_arset(p,&leadlag);
        //printf("p %s leadlag %d\n",p,leadlag);
       for (l=0; l<fdim; l++) {
             if(varindex==2) {
              ha_calvar[ha_calvarsize].Var2ADims[l]=0;
              ha_calvar[ha_calvarsize].Var2leadlag[l]=0;
              ha_calvar[ha_calvarsize].Var2SupSet[l]=0;
            } else {
              ha_calvar[ha_calvarsize].Var1ADims[l]=0;
              ha_calvar[ha_calvarsize].Var1leadlag[l]=0;
              ha_calvar[ha_calvarsize].Var1SupSet[l]=0;
            }
          if (strcmp(p,arSet[l].arIndx)==0) {
            if(varindex==2) {
              //ha_calvar[ha_calvarsize].Var2Type=1;
              //ha_calvar[ha_calvarsize].Var2BegAdd=ha_cof[index].begadd;
              //ha_calvar[*ha_calvarsize].Var2Size=fdim;//1;
              if (ha_set[ha_cof[index].setid[0]].size>ha_set[arSet[l].setid].size) {
                ha_calvar[ha_calvarsize].Var2SupSet[l]=1;
                for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[0]){ha_calvar[ha_calvarsize].Var2SSIndx[l]=sup;break;}
              }
              //ha_calvar[*ha_calvarsize].Var2SetBAdd[l]=arSet[l].SetBegAdd;
              ha_calvar[ha_calvarsize].Var2ADims[l]=ha_cof[index].antidims[0];
              ha_calvar[ha_calvarsize].Var2leadlag[l]=leadlag;
            } else {
              //ha_calvar[ha_calvarsize].Var1Type=1;
              //ha_calvar[ha_calvarsize].Var1BegAdd=ha_cof[index].begadd;
              //ha_calvar[*ha_calvarsize].Var1Size=fdim;//1;
              if (ha_set[ha_cof[index].setid[0]].size>ha_set[arSet[l].setid].size) {
                ha_calvar[ha_calvarsize].Var1SupSet[l]=1;
                for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[0]){ha_calvar[ha_calvarsize].Var1SSIndx[l]=sup;break;}
              }
              //ha_calvar[*ha_calvarsize].Var1SetBAdd[l]=arSet[l].SetBegAdd;
              ha_calvar[ha_calvarsize].Var1ADims[l]=ha_cof[index].antidims[0];
              ha_calvar[ha_calvarsize].Var1leadlag[l]=leadlag;
            }
          }
        }
        break;
      default:
        //i2=0;
        for (l1=0; l1<fdim; l1++) {
          if(varindex==2) {
            ha_calvar[ha_calvarsize].Var2ADims[l1]=0;
            ha_calvar[ha_calvarsize].Var2leadlag[l1]=0;
            ha_calvar[ha_calvarsize].Var2SupSet[l1]=0;
          } else {
            ha_calvar[ha_calvarsize].Var1ADims[l1]=0;
            ha_calvar[ha_calvarsize].Var1leadlag[l1]=0;
            ha_calvar[ha_calvarsize].Var1SupSet[l1]=0;
          }
        }
        for (l=0; l<ha_cof[index].size-1; l++) {
          p=strtok(NULL,",");
          leadlag=0;
          hnew_arset(p,&leadlag);
          for (l1=0; l1<fdim; l1++) {
            if (strcmp(p,arSet[l1].arIndx)==0) {
              //printf("p %s arset %s l %d l1 %d begadd %d arset %d\n",p,arSet[l1].arIndx,l,l1,ha_cof[index].begadd,arSet[l1].setid);
              if(varindex==2) {
                //ha_calvar[ha_calvarsize].Var2Type=1;
                //ha_calvar[ha_calvarsize].Var2BegAdd=ha_cof[index].begadd;
                //ha_calvar[*ha_calvarsize].Var2Size=fdim;//ha_cof[index].size;
                if (ha_set[ha_cof[index].setid[l]].size>ha_set[arSet[l1].setid].size) {
                  ha_calvar[ha_calvarsize].Var2SupSet[l1]=1;
                  for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l1].setid].subsetid[sup]==ha_cof[index].setid[l]){ha_calvar[ha_calvarsize].Var2SSIndx[l1]=sup;break;}
                }
                //ha_calvar[*ha_calvarsize].Var2SetBAdd[l1]=arSet[l1].SetBegAdd;
                ha_calvar[ha_calvarsize].Var2ADims[l1]=ha_cof[index].antidims[l];
                //printf("anti dim %d ",ha_calvar[ha_calvarsize].Var2ADims[l1]);
                ha_calvar[ha_calvarsize].Var2leadlag[l1]=leadlag;
              } else {
                //ha_calvar[ha_calvarsize].Var1Type=1;
                //ha_calvar[ha_calvarsize].Var1BegAdd=ha_cof[index].begadd;
                //ha_calvar[*ha_calvarsize].Var1Size=fdim;//ha_cof[index].size;
                if (ha_set[ha_cof[index].setid[l]].size>ha_set[arSet[l1].setid].size) {
                  ha_calvar[ha_calvarsize].Var1SupSet[l1]=1;
                  for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l1].setid].subsetid[sup]==ha_cof[index].setid[l]){ha_calvar[ha_calvarsize].Var1SSIndx[l1]=sup;break;}
                }
                //ha_calvar[*ha_calvarsize].Var1SetBAdd[l1]=arSet[l1].SetBegAdd;
                ha_calvar[ha_calvarsize].Var1ADims[l1]=ha_cof[index].antidims[l];
                //printf("anti dim %d ",ha_calvar[ha_calvarsize].Var1ADims[l1]);
                ha_calvar[ha_calvarsize].Var1leadlag[l1]=leadlag;
              }
              break;
            }
          }
        }
        p=strtok(NULL,"}");
        leadlag=0;
        hnew_arset(p,&leadlag);
        for (l1=0; l1<fdim; l1++) {
          if (strcmp(p,arSet[l1].arIndx)==0) {
            //printf("p %s arset %s l %d l1 %d dim %d arset %d\n",p,arSet[l1].arIndx,l,l1,ha_cof[index].dims[l],arSet[l1].SetSize);
            //if (ha_cof[index].dims[l]>arSet[l1].SetSize) i1=ha_setele[arSet[l1].SetBegAdd+arSet[l1].indx].setsh;
            //else i1=arSet[l1].indx;
            if(varindex==2) {
              if (ha_set[ha_cof[index].setid[l]].size>ha_set[arSet[l1].setid].size) {
                ha_calvar[ha_calvarsize].Var2SupSet[l1]=1;
                for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l1].setid].subsetid[sup]==ha_cof[index].setid[l]){ha_calvar[ha_calvarsize].Var2SSIndx[l1]=sup;break;}
              }
              //ha_calvar[*ha_calvarsize].Var2SetBAdd[l1]=arSet[l1].SetBegAdd;
              ha_calvar[ha_calvarsize].Var2ADims[l1]=ha_cof[index].antidims[l];
              //printf("anti dim %d\n",ha_calvar[ha_calvarsize].Var2ADims[l1]);
              ha_calvar[ha_calvarsize].Var2leadlag[l1]=leadlag;
            } else {
              if (ha_set[ha_cof[index].setid[l]].size>ha_set[arSet[l1].setid].size) {
                ha_calvar[ha_calvarsize].Var1SupSet[l1]=1;
                for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l1].setid].subsetid[sup]==ha_cof[index].setid[l]){ha_calvar[ha_calvarsize].Var1SSIndx[l1]=sup;break;}
              }
              //ha_calvar[*ha_calvarsize].Var1SetBAdd[l1]=arSet[l1].SetBegAdd;
              ha_calvar[ha_calvarsize].Var1ADims[l1]=ha_cof[index].antidims[l];
              //printf("anti dim %d\n",ha_calvar[ha_calvarsize].Var1ADims[l1]);
              ha_calvar[ha_calvarsize].Var1leadlag[l1]=leadlag;
            }
            break;
          }
        }
        break;
      }
      return 1;
    }
  } while (index--);

  index=nvar-1;
  do {
    if (strcmp(ha_var[index].cofname,p)==0) {
      //if(!ha_var[index].suplval)printf("Warning!!!! coefficient %s has not been supplied with values!\n",ha_var[index].cofname);
      //printf("p10 %s varname %s begadd %d ncof %d\n",p,ha_var[index].cofname,ha_var[index].begadd,ncofele);
      if(varindex==2) {
        ha_calvar[ha_calvarsize].Var2BegAdd=ha_var[index].begadd+ncofele;
        if(IsChange) ha_calvar[ha_calvarsize].Var2Type=6;
        else ha_calvar[ha_calvarsize].Var2Type=0;
      } else {
        ha_calvar[ha_calvarsize].Var1BegAdd=ha_var[index].begadd+ncofele;
        if(IsChange) ha_calvar[ha_calvarsize].Var1Type=6;
        else ha_calvar[ha_calvarsize].Var1Type=0;
      }
      //printf("varbegadd %d\n",ha_var[index].begadd);
      switch(ha_var[index].size) {
      case 0:
        for (l=0; l<fdim; l++) {
          if(varindex==2) {
            ha_calvar[ha_calvarsize].Var2ADims[l]=0;
            ha_calvar[ha_calvarsize].Var2leadlag[l]=0;
          } else {
            ha_calvar[ha_calvarsize].Var1ADims[l]=0;
            ha_calvar[ha_calvarsize].Var1leadlag[l]=0;
          }
        }
        break;
      case 1:
        p=strtok(NULL,"}");
        leadlag=0;
        hnew_arset(p,&leadlag);
        for (l=0; l<fdim; l++) {
             if(varindex==2) {
              ha_calvar[ha_calvarsize].Var2ADims[l]=0;
              ha_calvar[ha_calvarsize].Var2leadlag[l]=0;
              ha_calvar[ha_calvarsize].Var2SupSet[l]=0;
            } else {
              ha_calvar[ha_calvarsize].Var1ADims[l]=0;
              ha_calvar[ha_calvarsize].Var1leadlag[l]=0;
              ha_calvar[ha_calvarsize].Var1SupSet[l]=0;
            }
          if (strcmp(p,arSet[l].arIndx)==0) {
            if(varindex==2) {
              //ha_calvar[ha_calvarsize].Var2Type=0;
              //ha_calvar[ha_calvarsize].Var2BegAdd=ha_var[index].begadd;
              //ha_calvar[*ha_calvarsize].Var2Size=fdim;//1;
              if (ha_set[ha_var[index].setid[0]].size>ha_set[arSet[l].setid].size) {
                ha_calvar[ha_calvarsize].Var2SupSet[l]=1;
                for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[0]){ha_calvar[ha_calvarsize].Var2SSIndx[l]=sup;break;}
              }
              //ha_calvar[*ha_calvarsize].Var2SetBAdd[l]=arSet[l].SetBegAdd;
              ha_calvar[ha_calvarsize].Var2ADims[l]=ha_var[index].antidims[0];
              ha_calvar[ha_calvarsize].Var2leadlag[l]=leadlag;
            } else {
              //ha_calvar[ha_calvarsize].Var1Type=0;
              //ha_calvar[ha_calvarsize].Var1BegAdd=ha_var[index].begadd;
              //ha_calvar[*ha_calvarsize].Var1Size=fdim;//1;
              if (ha_set[ha_var[index].setid[0]].size>ha_set[arSet[l].setid].size) {
                ha_calvar[ha_calvarsize].Var1SupSet[l]=1;
                for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[0]){ha_calvar[ha_calvarsize].Var1SSIndx[l]=sup;break;}
              }
              //ha_calvar[*ha_calvarsize].Var1SetBAdd[l]=arSet[l].SetBegAdd;
              ha_calvar[ha_calvarsize].Var1ADims[l]=ha_var[index].antidims[0];
              ha_calvar[ha_calvarsize].Var1leadlag[l]=leadlag;
            }
          }
        }
        break;
      default:
        //i2=0;
        for (l1=0; l1<fdim; l1++) {
          if(varindex==2) {
            ha_calvar[ha_calvarsize].Var2ADims[l1]=0;
            ha_calvar[ha_calvarsize].Var2leadlag[l1]=0;
            ha_calvar[ha_calvarsize].Var2SupSet[l1]=0;
          } else {
            ha_calvar[ha_calvarsize].Var1ADims[l1]=0;
            ha_calvar[ha_calvarsize].Var1leadlag[l1]=0;
            ha_calvar[ha_calvarsize].Var1SupSet[l1]=0;
          }
        }
        for (l=0; l<ha_var[index].size-1; l++) {
          p=strtok(NULL,",");
          leadlag=0;
          hnew_arset(p,&leadlag);
          for (l1=0; l1<fdim; l1++) {
            if (strcmp(p,arSet[l1].arIndx)==0) {
              if(varindex==2) {
                //ha_calvar[ha_calvarsize].Var2Type=0;
                //ha_calvar[ha_calvarsize].Var2BegAdd=ha_var[index].begadd;
                //ha_calvar[*ha_calvarsize].Var2Size=fdim;//ha_var[index].size;
                if (ha_set[ha_var[index].setid[l]].size>ha_set[arSet[l1].setid].size) {
                  ha_calvar[ha_calvarsize].Var2SupSet[l1]=1;
                  for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l1].setid].subsetid[sup]==ha_var[index].setid[l]){ha_calvar[ha_calvarsize].Var2SSIndx[l1]=sup;break;}
                }
                //ha_calvar[*ha_calvarsize].Var2SetBAdd[l1]=arSet[l1].SetBegAdd;
                ha_calvar[ha_calvarsize].Var2ADims[l1]=ha_var[index].antidims[l];
                ha_calvar[ha_calvarsize].Var2leadlag[l1]=leadlag;
              } else {
                //ha_calvar[ha_calvarsize].Var1Type=0;
                //ha_calvar[ha_calvarsize].Var1BegAdd=ha_var[index].begadd;
                //ha_calvar[*ha_calvarsize].Var1Size=fdim;
                if (ha_set[ha_var[index].setid[l]].size>ha_set[arSet[l1].setid].size) {
                  ha_calvar[ha_calvarsize].Var1SupSet[l1]=1;
                  for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l1].setid].subsetid[sup]==ha_var[index].setid[l]){ha_calvar[ha_calvarsize].Var1SSIndx[l1]=sup;break;}
                }
                //ha_calvar[*ha_calvarsize].Var1SetBAdd[l1]=arSet[l1].SetBegAdd;
                ha_calvar[ha_calvarsize].Var1ADims[l1]=ha_var[index].antidims[l];
                ha_calvar[ha_calvarsize].Var1leadlag[l1]=leadlag;
              }
              break;
            }
          }
        }
        p=strtok(NULL,"}");
        leadlag=0;
        hnew_arset(p,&leadlag);
        for (l1=0; l1<fdim; l1++) {
          if (strcmp(p,arSet[l1].arIndx)==0) {
            //if (ha_var[index].dims[l]>arSet[l1].SetSize) i1=ha_setele[arSet[l1].SetBegAdd+arSet[l1].indx].setsh;
            //else i1=arSet[l1].indx;
            //i2=i2+i1*ha_var[index].antidims[l];
            if(varindex==2) {
              if (ha_set[ha_var[index].setid[l]].size>ha_set[arSet[l1].setid].size) {
                ha_calvar[ha_calvarsize].Var2SupSet[l1]=1;
                for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l1].setid].subsetid[sup]==ha_var[index].setid[l]){ha_calvar[ha_calvarsize].Var2SSIndx[l1]=sup;break;}
              }
              //ha_calvar[*ha_calvarsize].Var2SetBAdd[l1]=arSet[l1].SetBegAdd;
              ha_calvar[ha_calvarsize].Var2ADims[l1]=ha_var[index].antidims[l];
              ha_calvar[ha_calvarsize].Var2leadlag[l1]=leadlag;
            } else {
              if (ha_set[ha_var[index].setid[l]].size>ha_set[arSet[l1].setid].size) {
                ha_calvar[ha_calvarsize].Var1SupSet[l1]=1;
                for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l1].setid].subsetid[sup]==ha_var[index].setid[l]){ha_calvar[ha_calvarsize].Var1SSIndx[l1]=sup;break;}
              }
              //ha_calvar[*ha_calvarsize].Var1SetBAdd[l1]=arSet[l1].SetBegAdd;
              ha_calvar[ha_calvarsize].Var1ADims[l1]=ha_var[index].antidims[l];
              ha_calvar[ha_calvarsize].Var1leadlag[l1]=leadlag;
            }
            break;
          }
        }
        break;
      }
      return 1;
    }
  } while (index--);

  for (index=totalsum-1; index>-1; index--) {
    if (strcmp(sum_cof[index].sumname,p)==0) {
      if(varindex==2) {
        ha_calvar[ha_calvarsize].Var2BegAdd=sum_cof[index].begadd;
        ha_calvar[ha_calvarsize].Var2Type=2;
        for (l1=0; l1<sum_cof[index].size; l1++)ha_calvar[ha_calvarsize].Var2SSIndx[l1]=0;
      } else {
        ha_calvar[ha_calvarsize].Var1BegAdd=sum_cof[index].begadd;
        ha_calvar[ha_calvarsize].Var1Type=2;
        for (l1=0; l1<sum_cof[index].size; l1++)ha_calvar[ha_calvarsize].Var1SSIndx[l1]=0;
      }
      //printf("indx %d sumcof %s begadd %d upvar1 %d upvar2 %d\n",index,p,sum_cof[index].begadd,ha_calvar[ha_calvarsize].Var1SSIndx[1],ha_calvar[ha_calvarsize].Var2SSIndx[2]);
      switch(sum_cof[index].size) {
      case 0:
        for (l1=0; l1<fdim; l1++) {
          if(varindex==2) ha_calvar[ha_calvarsize].Var2ADims[l1]=0;
          else ha_calvar[ha_calvarsize].Var1ADims[l1]=0;
        }
        break;
      case 1:
        p=strtok(NULL,"}");
        for (l=0; l<fdim; l++) {
          if (strcmp(p,arSet[l].arIndx)==0) {
            if(varindex==2) {
              ha_calvar[ha_calvarsize].Var2ADims[l]=sum_cof[index].antidims[0];
            } else {
              ha_calvar[ha_calvarsize].Var1ADims[l]=sum_cof[index].antidims[0];
            }
          } else {
            if(varindex==2) {
              ha_calvar[ha_calvarsize].Var2ADims[l]=0;
            } else {
              ha_calvar[ha_calvarsize].Var1ADims[l]=0;
            }
          }
        }
        break;
      default:
        //i2=0;
        for (l1=0; l1<fdim; l1++) {
          if(varindex==2) ha_calvar[ha_calvarsize].Var2ADims[l1]=0;
          else ha_calvar[ha_calvarsize].Var1ADims[l1]=0;
        }
        for (l=0; l<sum_cof[index].size-1; l++) {
          p=strtok(NULL,",");
          for (l1=0; l1<fdim; l1++) {
            if (strcmp(p,arSet[l1].arIndx)==0) {
              //printf("p %s arset %s l %d l1 %d begadd %d arset %d\n",p,arSet[l1].arIndx,l,l1,sum_cof[index].size,arSet[l1].setid);
              //i1=arSet[l1].indx;
              //i2=i2+i1*sum_cof[index].antidims[l];
              if(varindex==2) {
                //ha_calvar[ha_calvarsize].Var2Type=2;
                //ha_calvar[ha_calvarsize].Var2BegAdd=sum_cof[index].begadd;
                //ha_calvar[*ha_calvarsize].Var2Size=fdim;//ha_var[index].size;
                //if (ha_var[index].dims[l]>arSet[l1].SetSize) ha_calvar[*ha_calvarsize].Var2SupSet[l1]=1;
                //ha_calvar[*ha_calvarsize].Var2SetBAdd[l1]=arSet[l1].SetBegAdd;
                ha_calvar[ha_calvarsize].Var2ADims[l1]=sum_cof[index].antidims[l];
              } else {
                //ha_calvar[ha_calvarsize].Var1Type=2;
                //ha_calvar[ha_calvarsize].Var1BegAdd=sum_cof[index].begadd;
                //ha_calvar[*ha_calvarsize].Var1Size=fdim;
                //if (ha_var[index].dims[l]>arSet[l1].SetSize) ha_calvar[*ha_calvarsize].Var1SupSet[l]=1;
                //ha_calvar[*ha_calvarsize].Var1SetBAdd[l1]=arSet[l1].SetBegAdd;
                ha_calvar[ha_calvarsize].Var1ADims[l1]=sum_cof[index].antidims[l];
              }
              break;
            }
          }
        }
        p=strtok(NULL,"}");
        for (l1=0; l1<fdim; l1++) {
          if (strcmp(p,arSet[l1].arIndx)==0) {
            //printf("p %s arset %s l %d l1 %d begadd %d arset %d\n",p,arSet[l1].arIndx,l,l1,sum_cof[index].size,arSet[l1].setid);
            //i1=arSet[l1].indx;
            //i2=i2+i1*sum_cof[index].antidims[l];
            if(varindex==2) {
              //if (ha_var[index].dims[l]>arSet[l1].SetSize) ha_calvar[*ha_calvarsize].Var2SupSet[l1]=1;
              //ha_calvar[*ha_calvarsize].Var2SetBAdd[l1]=arSet[l1].SetBegAdd;
              ha_calvar[ha_calvarsize].Var2ADims[l1]=sum_cof[index].antidims[l];
            } else {
              //if (ha_var[index].dims[l]>arSet[l1].SetSize) ha_calvar[*ha_calvarsize].Var1SupSet[l]=1;
              //ha_calvar[*ha_calvarsize].Var1SetBAdd[l1]=arSet[l1].SetBegAdd;
              ha_calvar[ha_calvarsize].Var1ADims[l1]=sum_cof[index].antidims[l];
            }
            break;
          }
        }
        break;
      }
      return 1;
    }
  }

  for (l1=0; l1<ha_calvarsize; l1++) {
    if (strcmp(var2,ha_calvar[l1].TmpVarName)==0) {
      if(varindex==2) {
        ha_calvar[ha_calvarsize].Var2Type=4;
        ha_calvar[ha_calvarsize].Var2BegAdd=l1;
      } else {
        ha_calvar[ha_calvarsize].Var1Type=4;
        ha_calvar[ha_calvarsize].Var1BegAdd=l1;
      }
      //break;
      return 1;
    }
  }
  //printf("var %s\n",var2);
  if (var2[0]>='0'&&var2[0]<='9') {
    //printf("var %s varindex %d\n",var2,varindex);
    if(varindex==2) {
      ha_calvar[ha_calvarsize].Var2Type=5;
      ha_calvar[ha_calvarsize].Var2Val=atof(var2);
    } else {
      ha_calvar[ha_calvarsize].Var1Type=5;
      ha_calvar[ha_calvarsize].Var1Val=atof(var2);
    }
    //printf("var %s varindex %d type %d val %f\n",var2,varindex,ha_calvar[ha_calvarsize].Var1Type,ha_calvar[ha_calvarsize].Var1Val);
    return 1;
  }
  return 0;
}

int hnew_intrpl(char *line) {
  char *p1,*p2,*p3;
  p1=strchr(line,'{');
  while (p1!=NULL) {
    p3=strchr(p1,'}');
    p2=strchr(p1,'+');
    //printf("p2 %s\n",p2);
    if(p2!=NULL) {
      while(p3-p2>0) {
        *p2='#';
        p2=strchr(p1,'+');
        if(p2==NULL)p2=p3;
      }
    }
    //printf("line %s\n",line);
    p2=strchr(p1,'-');
    if(p2!=NULL) {
      while(p3-p2>0) {
        *p2='!';
        p2=strchr(p1,'-');
        if(p2==NULL)p2=p3;
      }
    }
    //printf("line3 %s\n",line);
    p1=strchr(p3,'{');
    //printf("p1 %s\n",p1);
  }
  return 1;
}

int hnew_arset(char *p,int *leadlag) {
  char *plussign,*minsign;
  //printf("p %s\n",p);
  plussign=strchr(p,'#');
  minsign=strchr(p,'!');
  if(plussign!=NULL) {
    *leadlag=atoi(plussign+1);
    *plussign='\0';
  }
  if(minsign!=NULL) {
    *leadlag=-atoi(minsign+1);
    *minsign='\0';
  }
  //printf("p %s plumin %d\n",p,*leadlag);
  return 1;
}

int ha_newfparse(char *fomulain, ha_cgeset *ha_set,hcge_cof *ha_cof, uvadd ncof, hcge_cof *ha_var,uvadd nvar,uvadd ncofele,hcge_sumcof *sum_cof,uvdim totalsum,hcge_calvars *ha_calvar,uvdim *ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim) {
  int npar=0,npow=0,nmul=0,ndiv=0,nplu=0,nmin=0,j;
  *ha_calvarsize=0;
  npar=ha_cgenchf(fomulain, ')');
  //printf("fpart2a %s fdim %d\n",fomulain,fdim);
  if (npar==0) {
    npow=ha_cgenchf(fomulain, '^');
    if (npow>0) {
      ha_newfppow(fomulain,ha_set,npow,0,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,ha_calvarsize,arSet,fdim);
    }
    //eval=ha_cgefppow(fomulain,record,sofrecord,vartemppar,npar,vartemppow,npow,ncofele,ha_cof,ncof,ha_var,nvar,ha_setele,sum_cof,totalsum,ha_sumele,nsumele,arSet,fdim);
    nmul=ha_cgenchf(fomulain, '*');
    ndiv=ha_cgenchf(fomulain, '/');
    nmul=nmul+ndiv;
    if (nmul>0) {
      ha_newfpmuldiv(fomulain,ha_set,nmul,0,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,ha_calvarsize,arSet,fdim);
    }
    nplu=ha_cgenchf(fomulain, '+');
    nmin=ha_cgenchf(fomulain, '-');
    nplu=nplu+nmin;
    if (nplu>0) {
      ha_newfpplumin(fomulain,ha_set,nplu,0,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,ha_calvarsize,arSet,fdim);
      //printf("for %s\n",fomulain);
    }
    //else {//enable if switch plu 01
    //printf("fpart2b %s fdim %d\n",fomulain,fdim);
      hnew_varrepl(fomulain,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,*ha_calvarsize,arSet,fdim,1);
      ha_calvar[*ha_calvarsize].Oper=0;
      ha_calvar[*ha_calvarsize].TmpVarName[0]='\0';
      //printf("fpart2 %s fdim %d\n",fomulain,fdim);
      *ha_calvarsize=*ha_calvarsize+1;
    //}
    return 1;
  }
  int i;
  char *p=NULL;
  char fpart1[TABREADLINE],fpart2[TABREADLINE],fpart3[TABREADLINE];
  char interchar[TABREADLINE],interchar1[TABREADLINE];

  //printf("formula %s\n",fomulain);
  for (i=1; i<npar+2; i++) {
    p=strchr(fomulain,')');
    if (p!=NULL) {
      strncpy(fpart1, fomulain, p-fomulain);
      fpart1[p-fomulain] = '\0';
      strcpy(fpart3,p+1);
      p=strrchr(fpart1,'(');
      if (p!=NULL) {
        strcpy(fpart2, p+1);
        strncpy(fpart1,fomulain,p-fpart1);
        fpart1[p-fpart1] = '\0';
      } else {
        printf("Error in Formula!\n");
        return 0;
      }
    } else {
      fpart1[0]='\0';
      strcpy(fpart2,fomulain);
      fpart3[0]='\0';
    }
    //printf("fpart2 %s\n",fpart2);
    //eval=ha_cgefsimpar(fpart2,record,sofrecord,vartemppar,npar,ncofele,ha_cof,ncof,ha_var,nvar,vartempsum,totalsum,vartemppow,vartempmuldiv,vartempplu,arSet,fdim);
    npow=ha_cgenchf(fpart2, '^');
    if (npow>0) {
      ha_newfppow(fpart2,ha_set,npow,i,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,ha_calvarsize,arSet,fdim);
    }
    nmul=ha_cgenchf(fpart2, '*');
    ndiv=ha_cgenchf(fpart2, '/');
    nmul=nmul+ndiv;
    //printf("for %s\n",fpart2);
    if (nmul>0) {
      ha_newfpmuldiv(fpart2,ha_set,nmul,i,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,ha_calvarsize,arSet,fdim);
    }
    //printf("fpart2a %s\n",fpart2);
    nplu=ha_cgenchf(fpart2, '+');
    nmin=ha_cgenchf(fpart2, '-');
    nplu=nplu+nmin;
    if (nplu>0) {
      ha_newfpplumin(fpart2,ha_set,nplu,i,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,ha_calvarsize,arSet,fdim);
      //printf("fpart2 %s\n",fpart2);
    }
    //else {//enable if switch plu 01
      //printf("fpart2as %s\n",fpart2);
      if(strpbrk(fpart2,"=<>")==NULL){
      hnew_varrepl(fpart2,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,*ha_calvarsize,arSet,fdim,1);
      ha_calvar[*ha_calvarsize].Oper=0;
      ha_calvar[*ha_calvarsize].TmpVarName[0]='\0';
      //ha_calvar[*ha_calvarsize].TmpVarType=1;
      *ha_calvarsize=*ha_calvarsize+1;
      }
    //}
    //printf("j %d oper %d tmpvar %s tmpval var1 %d var2 %d\n",*ha_calvarsize-1,ha_calvar[*ha_calvarsize-1].Oper,ha_calvar[*ha_calvarsize-1].TmpVarName,ha_calvar[*ha_calvarsize-1].TmpVarVal,ha_calvar[*ha_calvarsize-1].Var1Type,ha_calvar[*ha_calvarsize-1].Var2Type);
    //printf("eval %lf\n",eval);
    sprintf(interchar1, "%d", i);
    interchar[0]='\0';
    if (i<10) {
      strcat(interchar,"ha_cgepar000");
    }
    if (9<i&&i<100) {
      strcat(interchar,"ha_cgepar00");
    }
    if (99<i&&i<1000) {
      strcat(interchar,"ha_cgepar0");
    }
    if (999<i&&i<10000) {
      strcat(interchar,"ha_cgepar");
    }
    if (i>10000) {
      strcat(interchar,"ha_cgepar");
      printf("Warning: Too many parenthesises\n");
    }
    strcat(interchar,interchar1);
    //strcpy(vartemppar[i-1].varname,interchar);
    strcpy(ha_calvar[*ha_calvarsize].TmpVarName,interchar);
    //ha_calvar[*ha_calvarsize].TmpVarType=1;
    ha_calvar[*ha_calvarsize].Oper=0;
    ha_calvar[*ha_calvarsize].Var1BegAdd=*ha_calvarsize-1;
    j=strlen(fpart1);
    //printf("fpart1s %s\n",fpart1);
    ha_calvar[*ha_calvarsize].Var1Type=4;
    if (j>3) if (fpart1[j-1]=='1'&&fpart1[j-2]=='0'&&tolower((int)fpart1[j-3])=='d'&&tolower((int)fpart1[j-4])=='i') {
        ha_calvar[*ha_calvarsize].Var1Type=41;
        fpart1[j-4]='\0';
      }
    if (j>3) if (fpart1[j-1]=='e'&&fpart1[j-2]=='g'&&fpart1[j-3]=='o'&&fpart1[j-4]=='l') {
        ha_calvar[*ha_calvarsize].Var1Type=43;
        //printf("form %s\npart1 %s\n",fomulain,fpart1);
        fpart1[j-4]='\0';
      }
    if (j==3) if (fpart1[j-1]=='s'&&fpart1[j-2]=='b'&&fpart1[j-3]=='a') {
        ha_calvar[*ha_calvarsize].Var1Type=42;
        fpart1[j-3]='\0';
      }
    if (j>3) if (fpart1[j-1]=='s'&&fpart1[j-2]=='b'&&fpart1[j-3]=='a') if(j==3||fpart1[j-4]==' '||fpart1[j-4]=='('||fpart1[j-4]=='+'||fpart1[j-4]=='-'||fpart1[j-4]=='*'||fpart1[j-4]=='/'||fpart1[j-4]=='^'||fpart1[j-4]==',') {
          ha_calvar[*ha_calvarsize].Var1Type=42;
          fpart1[j-3]='\0';
        }
    if (j==2) if (fpart1[j-1]=='f'&&fpart1[j-2]=='i') {
        ha_newfpif(fpart2,ha_set,2,i,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,ha_calvarsize,arSet,fdim);
        fpart1[j-2]='\0';
      }
    if (j>2) if (fpart1[j-1]=='f'&&fpart1[j-2]=='i') if(j==2||fpart1[j-3]==' '||fpart1[j-3]=='('||fpart1[j-3]=='+'||fpart1[j-3]=='-'||fpart1[j-3]=='*'||fpart1[j-3]=='/'||fpart1[j-3]=='^'||fpart1[j-3]==',') {
          ha_newfpif(fpart2,ha_set,2,i,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,ha_calvarsize,arSet,fdim);
          fpart1[j-2]='\0';
        }
    //printf("var1type %d\n",ha_calvar[*ha_calvarsize].Var1Type);
    *ha_calvarsize=*ha_calvarsize+1;
    //printf("var1type %d\n",ha_calvar[*ha_calvarsize-1].Var1Type);
    //printf("cs %d\n",*ha_calvarsize);
    strcat(fpart1, interchar);
    strcat(fpart1, fpart3);
    strcpy(fomulain,fpart1);
    //printf("j %d oper %d tmpvar %s tmpval var1 %d var2 %d\n",*ha_calvarsize-1,ha_calvar[*ha_calvarsize-1].Oper,ha_calvar[*ha_calvarsize-1].TmpVarName,ha_calvar[*ha_calvarsize-1].Var1Type,ha_calvar[*ha_calvarsize-1].Var2Type);
  }
  return 1;
}

ha_cgetype ha_newfpcal(ha_cgevar *record,ha_cgeset *ha_set,ha_cgesetele *ha_setele,ha_cgesumele *ha_sumele,hcge_calvars *ha_calvar,int ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim, ha_cgetype zerodivide) {
  int i;
  uvdim j;
  uvadd l=0,l1=0;
  ha_cgetype eval1=0,eval2=0,eval3=0;
  for (i=0; i<ha_calvarsize; i++) {
    //printf("i %d\n",i);
    switch(ha_calvar[i].Oper) {
    case 0:
      if (ha_calvar[i].Var1Type==0) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
          //printf("leadlag %d index %d l %d\n",ha_calvar[i].Var1leadlag[j],arSet[j].indx,l);
          //printf("var1 sub %d setsh %d\n",ha_calvar[i].Var1SupSet[j],ha_setele[arSet[j].SetBegAdd+arSet[j].indx].setsh);
        }
        ha_calvar[i].TmpVarVal=record[ha_calvar[i].Var1BegAdd+l].varval;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        break;
      }
      if (ha_calvar[i].Var1Type==2) {
        l=0;
        for (j=0; j<fdim; j++) {
          //if(ha_calvar[i].Var1SupSet[j]==1) {
          //l+=ha_calvar[i].Var1ADims[j]*ha_setele[arSet[j].SetBegAdd+arSet[j].indx].setsh;
          //} else {
          l+=ha_calvar[i].Var1ADims[j]*arSet[j].indx;
          //}
        }
        ha_calvar[i].TmpVarVal=ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
        //printf("add %d sume %f\n",ha_calvar[i].Var1BegAdd+l,ha_sumele[ha_calvar[i].Var1BegAdd+l].varval);
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        break;
      }
      if (ha_calvar[i].Var1Type==4) {
        ha_calvar[i].TmpVarVal=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        break;
      }
      if (ha_calvar[i].Var1Type==41) {
        if(ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal==0)ha_calvar[i].TmpVarVal=1;
        else ha_calvar[i].TmpVarVal=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        break;
      }
      if (ha_calvar[i].Var1Type==42) {
        if(ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal<0)ha_calvar[i].TmpVarVal=-ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
        else ha_calvar[i].TmpVarVal=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        break;
      }
      if (ha_calvar[i].Var1Type==43) {
        //printf("var %f varlog %f i %d\n",ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal,log(ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal),i);
        ha_calvar[i].TmpVarVal=log(ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal);
        //printf("var %lf log %lf\n",ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal,ha_calvar[i].TmpVarVal);
        break;
      }
      if (ha_calvar[i].Var1Type==5) {
        ha_calvar[i].TmpVarVal=ha_calvar[i].Var1Val;
        break;
      }
      if (ha_calvar[i].Var1Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
          //printf("leadlag %d index %d l %d\n",ha_calvar[i].Var1leadlag[j],arSet[j].indx,l);
          //printf("var1 sub %d setsh %d\n",ha_calvar[i].Var1SupSet[j],ha_setele[arSet[j].SetBegAdd+arSet[j].indx].setsh);
        }
        ha_calvar[i].TmpVarVal=record[ha_calvar[i].Var1BegAdd+l].csolpupd;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        break;
      }
      break;
    case 1:
      if(ha_calvar[i].Var1Type<3||ha_calvar[i].Var1Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
        }
      }
      if(ha_calvar[i].Var2Type<3||ha_calvar[i].Var2Type==6) {
        l1=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var2SupSet[j]==1) {
            l1+=ha_calvar[i].Var2ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]]+ha_calvar[i].Var2leadlag[j]);
          } else {
            l1+=ha_calvar[i].Var2ADims[j]*(arSet[j].indx+ha_calvar[i].Var2leadlag[j]);
          }
        }
      }
      //printf("l %d varadd %d vartype %d\n",l,ha_calvar[i].Var1BegAdd,ha_calvar[i].Var1Type);
      if(ha_calvar[i].Var1Type==0) eval1=record[ha_calvar[i].Var1BegAdd+l].varval;
      if(ha_calvar[i].Var1Type==2) eval1=ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
      if(ha_calvar[i].Var1Type==4) eval1=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
      if(ha_calvar[i].Var1Type==5) eval1=ha_calvar[i].Var1Val;
      if(ha_calvar[i].Var1Type==6) eval1=record[ha_calvar[i].Var1BegAdd+l].csolpupd;
      if(ha_calvar[i].Var2Type==0) eval2=record[ha_calvar[i].Var2BegAdd+l1].varval;
      if(ha_calvar[i].Var2Type==2) eval2=ha_sumele[ha_calvar[i].Var2BegAdd+l1].varval;
      if(ha_calvar[i].Var2Type==4) eval2=ha_calvar[ha_calvar[i].Var2BegAdd].TmpVarVal;
      if(ha_calvar[i].Var2Type==5) eval2=ha_calvar[i].Var2Val;
      if(ha_calvar[i].Var2Type==6) eval2=record[ha_calvar[i].Var2BegAdd+l1].csolpupd;
      ha_calvar[i].TmpVarVal=eval1*eval2;
      //printf("var %f i %d eval1 %f eval2 %f\n",ha_calvar[i].TmpVarVal,i,eval1,eval2);
      break;
    case 2:
      if(ha_calvar[i].Var1Type<3||ha_calvar[i].Var1Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
        }
      }
      if(ha_calvar[i].Var2Type<3||ha_calvar[i].Var2Type==6) {
        l1=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var2SupSet[j]==1) {
            l1+=ha_calvar[i].Var2ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]]+ha_calvar[i].Var2leadlag[j]);
          } else {
            l1+=ha_calvar[i].Var2ADims[j]*(arSet[j].indx+ha_calvar[i].Var2leadlag[j]);
          }
        }
      }
      if(ha_calvar[i].Var1Type==0) eval1=record[ha_calvar[i].Var1BegAdd+l].varval;
      if(ha_calvar[i].Var1Type==2) eval1=ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
      if(ha_calvar[i].Var1Type==4) eval1=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
      if(ha_calvar[i].Var1Type==5) eval1=ha_calvar[i].Var1Val;
      if(ha_calvar[i].Var1Type==6) eval1=record[ha_calvar[i].Var1BegAdd+l].csolpupd;
      if(ha_calvar[i].Var2Type==0) eval2=record[ha_calvar[i].Var2BegAdd+l1].varval;
      if(ha_calvar[i].Var2Type==2) eval2=ha_sumele[ha_calvar[i].Var2BegAdd+l1].varval;
      if(ha_calvar[i].Var2Type==4) eval2=ha_calvar[ha_calvar[i].Var2BegAdd].TmpVarVal;
      if(ha_calvar[i].Var2Type==5) eval2=ha_calvar[i].Var2Val;
      if(ha_calvar[i].Var2Type==6) eval2=record[ha_calvar[i].Var2BegAdd+l1].csolpupd;
      if(eval2==0) {
        ha_calvar[i].TmpVarVal=zerodivide;
      } else {
        ha_calvar[i].TmpVarVal=eval1/eval2;
      }
      //printf("eval1 %lf eval2 %lf var %f i %d\n",eval1,eval2,ha_calvar[i].TmpVarVal,i);
      break;
    case 3:
      if(ha_calvar[i].Var1Type==0) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
        }
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+record[ha_calvar[i].Var1BegAdd+l].varval;
        eval1=record[ha_calvar[i].Var1BegAdd+l].varval;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        //break;
      }
      if(ha_calvar[i].Var1Type==2) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]];
            //printf("la %d setele0 %d setele1 %d\n",l-l1,ha_calvar[i].Var1SSIndx[0],ha_calvar[i].Var1SSIndx[1]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*arSet[j].indx;
            //printf("lb %d\n",l-l1);
          }
        }
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
        eval1=ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
        //printf("eval1 %f i %d\n",eval1,l);
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        //break;
      }
      if (ha_calvar[i].Var1Type==4) {
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
        eval1=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
        //break;
      }
      if (ha_calvar[i].Var1Type==5) {
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+ha_calvar[i].Var1Val;
        eval1=ha_calvar[i].Var1Val;
        //break;
      }
      if(ha_calvar[i].Var1Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
        }
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+record[ha_calvar[i].Var1BegAdd+l].csolpupd;
        eval1=record[ha_calvar[i].Var1BegAdd+l].csolpupd;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        //break;
      }
      //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
      if(ha_calvar[i].Var2Type==0) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var2SupSet[j]==1) {
            l+=ha_calvar[i].Var2ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]]+ha_calvar[i].Var2leadlag[j]);
          } else {
            l+=ha_calvar[i].Var2ADims[j]*(arSet[j].indx+ha_calvar[i].Var2leadlag[j]);
          }
        }
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+record[ha_calvar[i].Var2BegAdd+l].varval;
        eval2=record[ha_calvar[i].Var2BegAdd+l].varval;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        //break;
      }
      if(ha_calvar[i].Var2Type==2) {
        l=0;
        for (j=0; j<fdim; j++) {
            l1=l;
          if(ha_calvar[i].Var2SupSet[j]==1) {
            l+=ha_calvar[i].Var2ADims[j]*ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]];
            //printf("j %d la %d setele0 %d setele1 %d setele2 %d\n",j,l-l1,ha_calvar[i].Var2SSIndx[0],ha_calvar[i].Var2SSIndx[1],ha_calvar[i].Var2SSIndx[2]);
          } else {
            l+=ha_calvar[i].Var2ADims[j]*arSet[j].indx;
            //printf("lb %d\n",l-l1);
          }
        }
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+ha_sumele[ha_calvar[i].Var2BegAdd+l].varval;
        eval2=ha_sumele[ha_calvar[i].Var2BegAdd+l].varval;
        //printf("eval2 %f i %d\n",eval2,l);
        //break;
      }
      if (ha_calvar[i].Var2Type==4) {
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+ha_calvar[ha_calvar[i].Var2BegAdd].TmpVarVal;
        eval2=ha_calvar[ha_calvar[i].Var2BegAdd].TmpVarVal;
        //break;
      }
      if (ha_calvar[i].Var2Type==5) {
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+ha_calvar[i].Var2Val;
        eval2=ha_calvar[i].Var2Val;
        //break;
      }
      if(ha_calvar[i].Var2Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var2SupSet[j]==1) {
            l+=ha_calvar[i].Var2ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]]+ha_calvar[i].Var2leadlag[j]);
          } else {
            l+=ha_calvar[i].Var2ADims[j]*(arSet[j].indx+ha_calvar[i].Var2leadlag[j]);
          }
        }
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+record[ha_calvar[i].Var2BegAdd+l].csolpupd;
        eval2=record[ha_calvar[i].Var2BegAdd+l].csolpupd;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        //break;
      }
      //if(ha_calvarsize==10)printf("var1 %lf var2 %lf\n",eval1,eval2);
      ha_calvar[i].TmpVarVal=eval1+eval2;
      break;
    case 4:
      if(ha_calvar[i].Var1Type==0) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
        }
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-record[ha_calvar[i].Var1BegAdd+l].varval;
        eval1=record[ha_calvar[i].Var1BegAdd+l].varval;
        //break;
      }
      if(ha_calvar[i].Var1Type==2) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]];
          } else {
            l+=ha_calvar[i].Var1ADims[j]*arSet[j].indx;
          }
        }
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
        eval1=ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
        //break;
      }
      if (ha_calvar[i].Var1Type==4) {
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
        eval1=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
        //break;
      }
      if (ha_calvar[i].Var1Type==5) {
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-ha_calvar[i].Var1Val;
        eval1=ha_calvar[i].Var1Val;
        //break;
      }
      if(ha_calvar[i].Var1Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
        }
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-record[ha_calvar[i].Var1BegAdd+l].csolpupd;
        eval1=record[ha_calvar[i].Var1BegAdd+l].csolpupd;
        //break;
      }
      //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
      if(ha_calvar[i].Var2Type==0) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var2SupSet[j]==1) {
            l+=ha_calvar[i].Var2ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]]+ha_calvar[i].Var2leadlag[j]);
          } else {
            l+=ha_calvar[i].Var2ADims[j]*(arSet[j].indx+ha_calvar[i].Var2leadlag[j]);
          }
        }
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-record[ha_calvar[i].Var2BegAdd+l].varval;
        eval2=record[ha_calvar[i].Var2BegAdd+l].varval;
        //break;
      }
      if(ha_calvar[i].Var2Type==2) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var2SupSet[j]==1) {
            l+=ha_calvar[i].Var2ADims[j]*ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]];
          } else {
            l+=ha_calvar[i].Var2ADims[j]*arSet[j].indx;
          }
        }
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-ha_sumele[ha_calvar[i].Var2BegAdd+l].varval;
        eval2=ha_sumele[ha_calvar[i].Var2BegAdd+l].varval;
        //break;
      }
      if (ha_calvar[i].Var2Type==4) {
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-ha_calvar[ha_calvar[i].Var2BegAdd].TmpVarVal;
        eval2=ha_calvar[ha_calvar[i].Var2BegAdd].TmpVarVal;
        //break;
      }
      if (ha_calvar[i].Var2Type==5) {
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-ha_calvar[i].Var2Val;
        eval2=ha_calvar[i].Var2Val;
        //break;
      }
      if(ha_calvar[i].Var2Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var2SupSet[j]==1) {
            l+=ha_calvar[i].Var2ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]]+ha_calvar[i].Var2leadlag[j]);
          } else {
            l+=ha_calvar[i].Var2ADims[j]*(arSet[j].indx+ha_calvar[i].Var2leadlag[j]);
          }
        }
        //ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-record[ha_calvar[i].Var2BegAdd+l].csolpupd;
        eval2=record[ha_calvar[i].Var2BegAdd+l].csolpupd;
        //break;
      }
      ha_calvar[i].TmpVarVal=eval1-eval2;
      break;
    case 5:
      if(ha_calvar[i].Var1Type<3||ha_calvar[i].Var1Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
        }
      }
      if(ha_calvar[i].Var2Type<3||ha_calvar[i].Var2Type==6) {
        l1=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var2SupSet[j]==1) {
            l1+=ha_calvar[i].Var2ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]]+ha_calvar[i].Var2leadlag[j]);
          } else {
            l1+=ha_calvar[i].Var2ADims[j]*(arSet[j].indx+ha_calvar[i].Var2leadlag[j]);
          }
        }
      }
      if(ha_calvar[i].Var1Type==0) eval1=record[ha_calvar[i].Var1BegAdd+l].varval;
      if(ha_calvar[i].Var1Type==2) eval1=ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
      if(ha_calvar[i].Var1Type==4) eval1=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
      if(ha_calvar[i].Var1Type==5) eval1=ha_calvar[i].Var1Val;
      if(ha_calvar[i].Var1Type==6) eval1=record[ha_calvar[i].Var1BegAdd+l].csolpupd;
      if(ha_calvar[i].Var2Type==0) eval2=record[ha_calvar[i].Var2BegAdd+l1].varval;
      if(ha_calvar[i].Var2Type==2) eval2=ha_sumele[ha_calvar[i].Var2BegAdd+l1].varval;
      if(ha_calvar[i].Var2Type==4) eval2=ha_calvar[ha_calvar[i].Var2BegAdd].TmpVarVal;
      if(ha_calvar[i].Var2Type==5) eval2=ha_calvar[i].Var2Val;
      if(ha_calvar[i].Var2Type==6) eval2=record[ha_calvar[i].Var2BegAdd+l1].csolpupd;
      if(eval1==0&&eval2<0) {
        ha_calvar[i].TmpVarVal=zerodivide;
      } else {
        if(eval1<0&&eval2-floor(eval2)!=0)printf("Serious errors: fraction power of negative number!!!!!!!");
        ha_calvar[i].TmpVarVal=pow(eval1,eval2);
      }
      //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
      break;
    default:
      if(ha_calvar[i].Var1Type==0) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
        }
        eval1=record[ha_calvar[i].Var1BegAdd+l].varval;
      }
      if(ha_calvar[i].Var1Type==2) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]];
          } else {
            l+=ha_calvar[i].Var1ADims[j]*arSet[j].indx;
          }
        }
        eval1=ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
      }
      if (ha_calvar[i].Var1Type==4) {
        eval1=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
      }
      if (ha_calvar[i].Var1Type==5) {
        eval1=ha_calvar[i].Var1Val;
      }
      if(ha_calvar[i].Var1Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
        }
        eval1=record[ha_calvar[i].Var1BegAdd+l].csolpupd;
      }
      if(ha_calvar[i].Var2Type==0) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var2SupSet[j]==1) {
            l+=ha_calvar[i].Var2ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]]+ha_calvar[i].Var2leadlag[j]);
          } else {
            l+=ha_calvar[i].Var2ADims[j]*(arSet[j].indx+ha_calvar[i].Var2leadlag[j]);
          }
        }
        eval2=record[ha_calvar[i].Var2BegAdd+l].varval;
      }
      if(ha_calvar[i].Var2Type==2) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var2SupSet[j]==1) {
            l+=ha_calvar[i].Var2ADims[j]*ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]];
          } else {
            l+=ha_calvar[i].Var2ADims[j]*arSet[j].indx;
          }
        }
        eval2=ha_sumele[ha_calvar[i].Var2BegAdd+l].varval;
      }
      if (ha_calvar[i].Var2Type==4) {
        eval2=ha_calvar[ha_calvar[i].Var2BegAdd].TmpVarVal;
      }
      if (ha_calvar[i].Var2Type==5) {
        eval2=ha_calvar[i].Var2Val;
      }
      if(ha_calvar[i].Var2Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var2SupSet[j]==1) {
            l+=ha_calvar[i].Var2ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]]+ha_calvar[i].Var2leadlag[j]);
          } else {
            l+=ha_calvar[i].Var2ADims[j]*(arSet[j].indx+ha_calvar[i].Var2leadlag[j]);
          }
        }
        eval2=record[ha_calvar[i].Var2BegAdd+l].csolpupd;
      }

      if(ha_calvar[i].Var3Type==0) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var3SupSet[j]==1) {
            l+=ha_calvar[i].Var3ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var3SSIndx[j]]+ha_calvar[i].Var3leadlag[j]);
          } else {
            l+=ha_calvar[i].Var3ADims[j]*(arSet[j].indx+ha_calvar[i].Var3leadlag[j]);
          }
        }
        eval3=record[ha_calvar[i].Var3BegAdd+l].varval;
      }
      if(ha_calvar[i].Var3Type==2) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var3SupSet[j]==1) {
            l+=ha_calvar[i].Var3ADims[j]*ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var3SSIndx[j]];
          } else {
            l+=ha_calvar[i].Var3ADims[j]*arSet[j].indx;
          }
        }
        eval3=ha_sumele[ha_calvar[i].Var3BegAdd+l].varval;
      }
      if (ha_calvar[i].Var3Type==4) {
        eval3=ha_calvar[ha_calvar[i].Var3BegAdd].TmpVarVal;
      }
      if (ha_calvar[i].Var3Type==5) {
        eval3=ha_calvar[i].Var3Val;
      }
      if(ha_calvar[i].Var3Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var3SupSet[j]==1) {
            l+=ha_calvar[i].Var3ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var3SSIndx[j]]+ha_calvar[i].Var3leadlag[j]);
          } else {
            l+=ha_calvar[i].Var3ADims[j]*(arSet[j].indx+ha_calvar[i].Var3leadlag[j]);
          }
        }
        eval3=record[ha_calvar[i].Var3BegAdd+l].csolpupd;
      }
      //printf("i %d ini %lf\n",i,ha_calvar[i].TmpVarVal);
      if(ha_calvar[i].Oper==71)if(eval1==eval2)ha_calvar[i].TmpVarVal=eval3;else ha_calvar[i].TmpVarVal=0;
      if(ha_calvar[i].Oper==72)if(eval1>eval2)ha_calvar[i].TmpVarVal=eval3;else ha_calvar[i].TmpVarVal=0;
      if(ha_calvar[i].Oper==73)if(eval1<eval2)ha_calvar[i].TmpVarVal=eval3;else ha_calvar[i].TmpVarVal=0;
      if(ha_calvar[i].Oper==74)if(eval1!=eval2)ha_calvar[i].TmpVarVal=eval3;else ha_calvar[i].TmpVarVal=0;
      if(ha_calvar[i].Oper==75)if(eval1<=eval2)ha_calvar[i].TmpVarVal=eval3;else ha_calvar[i].TmpVarVal=0;
      if(ha_calvar[i].Oper==76)if(eval1>=eval2)ha_calvar[i].TmpVarVal=eval3;else ha_calvar[i].TmpVarVal=0;
      break;
    }
    //printf("i %d var1 %d var2 %d val %lf\n",i,ha_calvar[i].Var1BegAdd,ha_calvar[i].Var2BegAdd,ha_calvar[i].TmpVarVal);
  }
  return ha_calvar[i-1].TmpVarVal;
}

ha_cgetype ha_newfpcal01(ha_cgevar *record,ha_cgeset *ha_set,ha_cgesetele *ha_setele,ha_cgesumele *ha_sumele,hcge_calvars *ha_calvar,int ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim, ha_cgetype zerodivide) {
  int i;
  uvdim j;
  uvadd l=0,l1=0;
  ha_cgetype eval1=0,eval2=0;
  for (i=0; i<ha_calvarsize; i++) {
    switch(ha_calvar[i].Oper) {
    case 0:
      if (ha_calvar[i].Var1Type==0) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
          //printf("leadlag %d index %d l %d\n",ha_calvar[i].Var1leadlag[j],arSet[j].indx,l);
          //printf("var1 sub %d setsh %d\n",ha_calvar[i].Var1SupSet[j],ha_setele[arSet[j].SetBegAdd+arSet[j].indx].setsh);
        }
        ha_calvar[i].TmpVarVal=record[ha_calvar[i].Var1BegAdd+l].varval;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        break;
      }
      if (ha_calvar[i].Var1Type==2) {
        l=0;
        for (j=0; j<fdim; j++) {
          //if(ha_calvar[i].Var1SupSet[j]==1) {
          //l+=ha_calvar[i].Var1ADims[j]*ha_setele[arSet[j].SetBegAdd+arSet[j].indx].setsh;
          //} else {
          l+=ha_calvar[i].Var1ADims[j]*arSet[j].indx;
          //}
        }
        ha_calvar[i].TmpVarVal=ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
        //printf("add %d sume %f\n",ha_calvar[i].Var1BegAdd+l,ha_sumele[ha_calvar[i].Var1BegAdd+l].varval);
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        break;
      }
      if (ha_calvar[i].Var1Type==4) {
        ha_calvar[i].TmpVarVal=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        break;
      }
      if (ha_calvar[i].Var1Type==41) {
        if(ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal==0)ha_calvar[i].TmpVarVal=1;
        else ha_calvar[i].TmpVarVal=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        break;
      }
      if (ha_calvar[i].Var1Type==42) {
        if(ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal<0)ha_calvar[i].TmpVarVal=-ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
        else ha_calvar[i].TmpVarVal=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        break;
      }
      if (ha_calvar[i].Var1Type==43) {
        //printf("var %f varlog %f i %d\n",ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal,log(ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal),i);
        ha_calvar[i].TmpVarVal=log(ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal);
        //printf("var %lf log %lf\n",ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal,ha_calvar[i].TmpVarVal);
        break;
      }
      if (ha_calvar[i].Var1Type==5) {
        ha_calvar[i].TmpVarVal=ha_calvar[i].Var1Val;
        break;
      }
      if (ha_calvar[i].Var1Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
          //printf("leadlag %d index %d l %d\n",ha_calvar[i].Var1leadlag[j],arSet[j].indx,l);
          //printf("var1 sub %d setsh %d\n",ha_calvar[i].Var1SupSet[j],ha_setele[arSet[j].SetBegAdd+arSet[j].indx].setsh);
        }
        ha_calvar[i].TmpVarVal=record[ha_calvar[i].Var1BegAdd+l].csolpupd;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        break;
      }
      break;
    case 1:
      if(ha_calvar[i].Var1Type<3||ha_calvar[i].Var1Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
        }
      }
      if(ha_calvar[i].Var2Type<3||ha_calvar[i].Var2Type==6) {
        l1=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var2SupSet[j]==1) {
            l1+=ha_calvar[i].Var2ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]]+ha_calvar[i].Var2leadlag[j]);
          } else {
            l1+=ha_calvar[i].Var2ADims[j]*(arSet[j].indx+ha_calvar[i].Var2leadlag[j]);
          }
        }
      }
      //printf("l %d varadd %d vartype %d\n",l,ha_calvar[i].Var1BegAdd,ha_calvar[i].Var1Type);
      if(ha_calvar[i].Var1Type==0) eval1=record[ha_calvar[i].Var1BegAdd+l].varval;
      if(ha_calvar[i].Var1Type==2) eval1=ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
      if(ha_calvar[i].Var1Type==4) eval1=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
      if(ha_calvar[i].Var1Type==5) eval1=ha_calvar[i].Var1Val;
      if(ha_calvar[i].Var1Type==6) eval1=record[ha_calvar[i].Var1BegAdd+l].csolpupd;
      if(ha_calvar[i].Var2Type==0) eval2=record[ha_calvar[i].Var2BegAdd+l1].varval;
      if(ha_calvar[i].Var2Type==2) eval2=ha_sumele[ha_calvar[i].Var2BegAdd+l1].varval;
      if(ha_calvar[i].Var2Type==4) eval2=ha_calvar[ha_calvar[i].Var2BegAdd].TmpVarVal;
      if(ha_calvar[i].Var2Type==5) eval2=ha_calvar[i].Var2Val;
      if(ha_calvar[i].Var2Type==6) eval2=record[ha_calvar[i].Var2BegAdd+l1].csolpupd;
      ha_calvar[i].TmpVarVal=eval1*eval2;
      //printf("var %f i %d eval1 %f eval2 %f\n",ha_calvar[i].TmpVarVal,i,eval1,eval2);
      break;
    case 2:
      if(ha_calvar[i].Var1Type<3||ha_calvar[i].Var1Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
        }
      }
      if(ha_calvar[i].Var2Type<3||ha_calvar[i].Var2Type==6) {
        l1=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var2SupSet[j]==1) {
            l1+=ha_calvar[i].Var2ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]]+ha_calvar[i].Var2leadlag[j]);
          } else {
            l1+=ha_calvar[i].Var2ADims[j]*(arSet[j].indx+ha_calvar[i].Var2leadlag[j]);
          }
        }
      }
      if(ha_calvar[i].Var1Type==0) eval1=record[ha_calvar[i].Var1BegAdd+l].varval;
      if(ha_calvar[i].Var1Type==2) eval1=ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
      if(ha_calvar[i].Var1Type==4) eval1=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
      if(ha_calvar[i].Var1Type==5) eval1=ha_calvar[i].Var1Val;
      if(ha_calvar[i].Var1Type==6) eval1=record[ha_calvar[i].Var1BegAdd+l].csolpupd;
      if(ha_calvar[i].Var2Type==0) eval2=record[ha_calvar[i].Var2BegAdd+l1].varval;
      if(ha_calvar[i].Var2Type==2) eval2=ha_sumele[ha_calvar[i].Var2BegAdd+l1].varval;
      if(ha_calvar[i].Var2Type==4) eval2=ha_calvar[ha_calvar[i].Var2BegAdd].TmpVarVal;
      if(ha_calvar[i].Var2Type==5) eval2=ha_calvar[i].Var2Val;
      if(ha_calvar[i].Var2Type==6) eval2=record[ha_calvar[i].Var2BegAdd+l1].csolpupd;
      if(eval2==0) {
        ha_calvar[i].TmpVarVal=zerodivide;
      } else {
        ha_calvar[i].TmpVarVal=eval1/eval2;
      }
      //printf("eval1 %lf eval2 %lf var %f i %d\n",eval1,eval2,ha_calvar[i].TmpVarVal,i);
      break;
    case 3:
      if(ha_calvar[i].Var1Type==0) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
        }
        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+record[ha_calvar[i].Var1BegAdd+l].varval;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        break;
      }
      if(ha_calvar[i].Var1Type==2) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]];
          } else {
            l+=ha_calvar[i].Var1ADims[j]*arSet[j].indx;
          }
        }
        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        break;
      }
      if (ha_calvar[i].Var1Type==4) {
        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
        break;
      }
      if (ha_calvar[i].Var1Type==5) {
        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+ha_calvar[i].Var1Val;
        break;
      }
      if(ha_calvar[i].Var1Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
        }
        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+record[ha_calvar[i].Var1BegAdd+l].csolpupd;
        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
        break;
      }
      //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
      break;
    case 4:
      if(ha_calvar[i].Var1Type==0) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
        }
        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-record[ha_calvar[i].Var1BegAdd+l].varval;
        break;
      }
      if(ha_calvar[i].Var1Type==2) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]];
          } else {
            l+=ha_calvar[i].Var1ADims[j]*arSet[j].indx;
          }
        }
        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
        break;
      }
      if (ha_calvar[i].Var1Type==4) {
        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
        break;
      }
      if (ha_calvar[i].Var1Type==5) {
        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-ha_calvar[i].Var1Val;
        break;
      }
      if(ha_calvar[i].Var1Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
        }
        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-record[ha_calvar[i].Var1BegAdd+l].csolpupd;
        break;
      }
      //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
      break;
    case 5:
      if(ha_calvar[i].Var1Type<3||ha_calvar[i].Var1Type==6) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var1SupSet[j]==1) {
            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
          } else {
            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
          }
        }
      }
      if(ha_calvar[i].Var2Type<3||ha_calvar[i].Var2Type==6) {
        l1=0;
        for (j=0; j<fdim; j++) {
          if(ha_calvar[i].Var2SupSet[j]==1) {
            l1+=ha_calvar[i].Var2ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]]+ha_calvar[i].Var2leadlag[j]);
          } else {
            l1+=ha_calvar[i].Var2ADims[j]*(arSet[j].indx+ha_calvar[i].Var2leadlag[j]);
          }
        }
      }
      if(ha_calvar[i].Var1Type==0) eval1=record[ha_calvar[i].Var1BegAdd+l].varval;
      if(ha_calvar[i].Var1Type==2) eval1=ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
      if(ha_calvar[i].Var1Type==4) eval1=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
      if(ha_calvar[i].Var1Type==5) eval1=ha_calvar[i].Var1Val;
      if(ha_calvar[i].Var1Type==6) eval1=record[ha_calvar[i].Var1BegAdd+l].csolpupd;
      if(ha_calvar[i].Var2Type==0) eval2=record[ha_calvar[i].Var2BegAdd+l1].varval;
      if(ha_calvar[i].Var2Type==2) eval2=ha_sumele[ha_calvar[i].Var2BegAdd+l1].varval;
      if(ha_calvar[i].Var2Type==4) eval2=ha_calvar[ha_calvar[i].Var2BegAdd].TmpVarVal;
      if(ha_calvar[i].Var2Type==5) eval2=ha_calvar[i].Var2Val;
      if(ha_calvar[i].Var2Type==6) eval2=record[ha_calvar[i].Var2BegAdd+l1].csolpupd;
      if(eval1==0&&eval2<0) {
        ha_calvar[i].TmpVarVal=zerodivide;
      } else {
        if(eval1<0&&eval2-floor(eval2)!=0)printf("Serious errors: fraction power of negative number!!!!!!!");
        ha_calvar[i].TmpVarVal=pow(eval1,eval2);
      }
      //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
      break;
    default:
      break;
    }
    //printf("res varval %f oper %d\n",ha_calvar[i].TmpVarVal,ha_calvar[i].Oper);
  }
  return ha_calvar[i-1].TmpVarVal;
}

//ha_cgetype ha_newfpcalshow(ha_cgevar *record,ha_cgeset *ha_set,ha_cgesetele *ha_setele,ha_cgesumele *ha_sumele,hcge_calvars *ha_calvar,int ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim, ha_cgetype zerodivide,int show) {
//  int i;
//  uvdim j;
//  uvadd l=0,l1=0;
//  ha_cgetype eval1=0,eval2=0;
//  for (i=0; i<ha_calvarsize; i++) {
//    switch(ha_calvar[i].Oper) {
//    case 0:
//      if (ha_calvar[i].Var1Type==0) {
//        l=0;
//        for (j=0; j<fdim; j++) {
//          if(ha_calvar[i].Var1SupSet[j]==1) {
//            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
//          } else {
//            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
//          }
//          //printf("leadlag %d index %d l %d\n",ha_calvar[i].Var1leadlag[j],arSet[j].indx,l);
//          //printf("var1 sub %d setsh %d\n",ha_calvar[i].Var1SupSet[j],ha_setele[arSet[j].SetBegAdd+arSet[j].indx].setsh);
//        }
//        ha_calvar[i].TmpVarVal=record[ha_calvar[i].Var1BegAdd+l].varval;
//        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
//        break;
//      }
//      if (ha_calvar[i].Var1Type==2) {
//        l=0;
//        for (j=0; j<fdim; j++) {
//          //if(ha_calvar[i].Var1SupSet[j]==1) {
//          //l+=ha_calvar[i].Var1ADims[j]*ha_setele[arSet[j].SetBegAdd+arSet[j].indx].setsh;
//          //} else {
//          l+=ha_calvar[i].Var1ADims[j]*arSet[j].indx;
//          //}
//        }
//        ha_calvar[i].TmpVarVal=ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
//        //printf("add %d sume %f\n",ha_calvar[i].Var1BegAdd+l,ha_sumele[ha_calvar[i].Var1BegAdd+l].varval);
//        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
//        break;
//      }
//      if (ha_calvar[i].Var1Type==4) {
//        ha_calvar[i].TmpVarVal=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
//        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
//        break;
//      }
//      if (ha_calvar[i].Var1Type==41) {
//        if(ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal==0)ha_calvar[i].TmpVarVal=1;
//        else ha_calvar[i].TmpVarVal=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
//        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
//        break;
//      }
//      if (ha_calvar[i].Var1Type==42) {
//        if(ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal<0)ha_calvar[i].TmpVarVal=-ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
//        else ha_calvar[i].TmpVarVal=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
//        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
//        break;
//      }
//      if (ha_calvar[i].Var1Type==43) {
//        //printf("var %f varlog %f i %d\n",ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal,log(ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal),i);
//        ha_calvar[i].TmpVarVal=log(ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal);
//        //printf("var %lf log %lf\n",ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal,ha_calvar[i].TmpVarVal);
//        break;
//      }
//      if (ha_calvar[i].Var1Type==5) {
//        ha_calvar[i].TmpVarVal=ha_calvar[i].Var1Val;
//        break;
//      }
//      if (ha_calvar[i].Var1Type==6) {
//        l=0;
//        for (j=0; j<fdim; j++) {
//          if(ha_calvar[i].Var1SupSet[j]==1) {
//            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
//          } else {
//            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
//          }
//          //printf("leadlag %d index %d l %d\n",ha_calvar[i].Var1leadlag[j],arSet[j].indx,l);
//          //printf("var1 sub %d setsh %d\n",ha_calvar[i].Var1SupSet[j],ha_setele[arSet[j].SetBegAdd+arSet[j].indx].setsh);
//        }
//        ha_calvar[i].TmpVarVal=record[ha_calvar[i].Var1BegAdd+l].csolpupd;
//        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
//        break;
//      }
//      break;
//    case 1:
//      if(ha_calvar[i].Var1Type<3||ha_calvar[i].Var1Type==6) {
//        l=0;
//        for (j=0; j<fdim; j++) {
//          if(ha_calvar[i].Var1SupSet[j]==1) {
//            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
//          } else {
//            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
//          }
//        }
//      }
//      if(ha_calvar[i].Var2Type<3||ha_calvar[i].Var2Type==6) {
//        l1=0;
//        for (j=0; j<fdim; j++) {
//          if(ha_calvar[i].Var2SupSet[j]==1) {
//            l1+=ha_calvar[i].Var2ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]]+ha_calvar[i].Var2leadlag[j]);
//          } else {
//            l1+=ha_calvar[i].Var2ADims[j]*(arSet[j].indx+ha_calvar[i].Var2leadlag[j]);
//          }
//        }
//      }
//      //printf("l %d varadd %d vartype %d\n",l,ha_calvar[i].Var1BegAdd,ha_calvar[i].Var1Type);
//      if(ha_calvar[i].Var1Type==0) eval1=record[ha_calvar[i].Var1BegAdd+l].varval;
//      if(ha_calvar[i].Var1Type==2) eval1=ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
//      if(ha_calvar[i].Var1Type==4) eval1=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
//      if(ha_calvar[i].Var1Type==5) eval1=ha_calvar[i].Var1Val;
//      if(ha_calvar[i].Var1Type==6) eval1=record[ha_calvar[i].Var1BegAdd+l].csolpupd;
//      if(ha_calvar[i].Var2Type==0) eval2=record[ha_calvar[i].Var2BegAdd+l1].varval;
//      if(ha_calvar[i].Var2Type==2) eval2=ha_sumele[ha_calvar[i].Var2BegAdd+l1].varval;
//      if(ha_calvar[i].Var2Type==4) eval2=ha_calvar[ha_calvar[i].Var2BegAdd].TmpVarVal;
//      if(ha_calvar[i].Var2Type==5) eval2=ha_calvar[i].Var2Val;
//      if(ha_calvar[i].Var2Type==6) eval2=record[ha_calvar[i].Var2BegAdd+l1].csolpupd;
//      ha_calvar[i].TmpVarVal=eval1*eval2;
//      printf("var %f i %d eval1 %f eval2 %f record %f record add %d \n",ha_calvar[i].TmpVarVal,i,eval1,eval2,record[ha_calvar[i].Var2BegAdd+l1].varval,ha_calvar[i].Var2BegAdd+l1);
//      break;
//    case 2:
//      if(ha_calvar[i].Var1Type<3||ha_calvar[i].Var1Type==6) {
//        l=0;
//        for (j=0; j<fdim; j++) {
//          if(ha_calvar[i].Var1SupSet[j]==1) {
//            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
//          } else {
//            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
//          }
//        }
//      }
//      if(ha_calvar[i].Var2Type<3||ha_calvar[i].Var2Type==6) {
//        l1=0;
//        for (j=0; j<fdim; j++) {
//          if(ha_calvar[i].Var2SupSet[j]==1) {
//            l1+=ha_calvar[i].Var2ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]]+ha_calvar[i].Var2leadlag[j]);
//          } else {
//            l1+=ha_calvar[i].Var2ADims[j]*(arSet[j].indx+ha_calvar[i].Var2leadlag[j]);
//          }
//        }
//      }
//      if(ha_calvar[i].Var1Type==0) eval1=record[ha_calvar[i].Var1BegAdd+l].varval;
//      if(ha_calvar[i].Var1Type==2) eval1=ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
//      if(ha_calvar[i].Var1Type==4) eval1=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
//      if(ha_calvar[i].Var1Type==5) eval1=ha_calvar[i].Var1Val;
//      if(ha_calvar[i].Var1Type==6) eval1=record[ha_calvar[i].Var1BegAdd+l].csolpupd;
//      if(ha_calvar[i].Var2Type==0) eval2=record[ha_calvar[i].Var2BegAdd+l1].varval;
//      if(ha_calvar[i].Var2Type==2) eval2=ha_sumele[ha_calvar[i].Var2BegAdd+l1].varval;
//      if(ha_calvar[i].Var2Type==4) eval2=ha_calvar[ha_calvar[i].Var2BegAdd].TmpVarVal;
//      if(ha_calvar[i].Var2Type==5) eval2=ha_calvar[i].Var2Val;
//      if(ha_calvar[i].Var2Type==6) eval2=record[ha_calvar[i].Var2BegAdd+l1].csolpupd;
//      if(eval2==0) {
//        ha_calvar[i].TmpVarVal=zerodivide;
//      } else {
//        ha_calvar[i].TmpVarVal=eval1/eval2;
//      }
//      //printf("eval1 %lf eval2 %lf var %f i %d\n",eval1,eval2,ha_calvar[i].TmpVarVal,i);
//      break;
//    case 3:
//      if(ha_calvar[i].Var1Type==0) {
//        l=0;
//        for (j=0; j<fdim; j++) {
//          if(ha_calvar[i].Var1SupSet[j]==1) {
//            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
//          } else {
//            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
//          }
//        }
//        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+record[ha_calvar[i].Var1BegAdd+l].varval;
//        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
//        break;
//      }
//      if(ha_calvar[i].Var1Type==2) {
//        l=0;
//        for (j=0; j<fdim; j++) {
//          if(ha_calvar[i].Var1SupSet[j]==1) {
//            l+=ha_calvar[i].Var1ADims[j]*ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]];
//          } else {
//            l+=ha_calvar[i].Var1ADims[j]*arSet[j].indx;
//          }
//        }
//        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
//        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
//        break;
//      }
//      if (ha_calvar[i].Var1Type==4) {
//        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
//        break;
//      }
//      if (ha_calvar[i].Var1Type==5) {
//        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+ha_calvar[i].Var1Val;
//        break;
//      }
//      if(ha_calvar[i].Var1Type==6) {
//        l=0;
//        for (j=0; j<fdim; j++) {
//          if(ha_calvar[i].Var1SupSet[j]==1) {
//            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
//          } else {
//            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
//          }
//        }
//        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal+record[ha_calvar[i].Var1BegAdd+l].csolpupd;
//        //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
//        break;
//      }
//      //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
//      break;
//    case 4:
//      if(ha_calvar[i].Var1Type==0) {
//        l=0;
//        for (j=0; j<fdim; j++) {
//          if(ha_calvar[i].Var1SupSet[j]==1) {
//            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
//          } else {
//            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
//          }
//        }
//        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-record[ha_calvar[i].Var1BegAdd+l].varval;
//        break;
//      }
//      if(ha_calvar[i].Var1Type==2) {
//        l=0;
//        for (j=0; j<fdim; j++) {
//          if(ha_calvar[i].Var1SupSet[j]==1) {
//            l+=ha_calvar[i].Var1ADims[j]*ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]];
//          } else {
//            l+=ha_calvar[i].Var1ADims[j]*arSet[j].indx;
//          }
//        }
//        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
//        break;
//      }
//      if (ha_calvar[i].Var1Type==4) {
//        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
//        break;
//      }
//      if (ha_calvar[i].Var1Type==5) {
//        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-ha_calvar[i].Var1Val;
//        break;
//      }
//      if(ha_calvar[i].Var1Type==6) {
//        l=0;
//        for (j=0; j<fdim; j++) {
//          if(ha_calvar[i].Var1SupSet[j]==1) {
//            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
//          } else {
//            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
//          }
//        }
//        ha_calvar[i].TmpVarVal=ha_calvar[i-1].TmpVarVal-record[ha_calvar[i].Var1BegAdd+l].csolpupd;
//        break;
//      }
//      //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
//      break;
//    case 5:
//      if(ha_calvar[i].Var1Type<3||ha_calvar[i].Var1Type==6) {
//        l=0;
//        for (j=0; j<fdim; j++) {
//          if(ha_calvar[i].Var1SupSet[j]==1) {
//            l+=ha_calvar[i].Var1ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var1SSIndx[j]]+ha_calvar[i].Var1leadlag[j]);
//          } else {
//            l+=ha_calvar[i].Var1ADims[j]*(arSet[j].indx+ha_calvar[i].Var1leadlag[j]);
//          }
//        }
//      }
//      if(ha_calvar[i].Var2Type<3||ha_calvar[i].Var2Type==6) {
//        l1=0;
//        for (j=0; j<fdim; j++) {
//          if(ha_calvar[i].Var2SupSet[j]==1) {
//            l1+=ha_calvar[i].Var2ADims[j]*(ha_setele[ha_set[arSet[j].setid].begadd+arSet[j].indx].setsh[ha_calvar[i].Var2SSIndx[j]]+ha_calvar[i].Var2leadlag[j]);
//          } else {
//            l1+=ha_calvar[i].Var2ADims[j]*(arSet[j].indx+ha_calvar[i].Var2leadlag[j]);
//          }
//        }
//      }
//      if(ha_calvar[i].Var1Type==0) eval1=record[ha_calvar[i].Var1BegAdd+l].varval;
//      if(ha_calvar[i].Var1Type==2) eval1=ha_sumele[ha_calvar[i].Var1BegAdd+l].varval;
//      if(ha_calvar[i].Var1Type==4) eval1=ha_calvar[ha_calvar[i].Var1BegAdd].TmpVarVal;
//      if(ha_calvar[i].Var1Type==5) eval1=ha_calvar[i].Var1Val;
//      if(ha_calvar[i].Var1Type==6) eval1=record[ha_calvar[i].Var1BegAdd+l].csolpupd;
//      if(ha_calvar[i].Var2Type==0) eval2=record[ha_calvar[i].Var2BegAdd+l1].varval;
//      if(ha_calvar[i].Var2Type==2) eval2=ha_sumele[ha_calvar[i].Var2BegAdd+l1].varval;
//      if(ha_calvar[i].Var2Type==4) eval2=ha_calvar[ha_calvar[i].Var2BegAdd].TmpVarVal;
//      if(ha_calvar[i].Var2Type==5) eval2=ha_calvar[i].Var2Val;
//      if(ha_calvar[i].Var2Type==6) eval2=record[ha_calvar[i].Var2BegAdd+l1].csolpupd;
//      if(eval1==0&&eval2<0) {
//        ha_calvar[i].TmpVarVal=zerodivide;
//      } else {
//        ha_calvar[i].TmpVarVal=pow(eval1,eval2);
//      }
//      //printf("var %f i %d\n",ha_calvar[i].TmpVarVal,i);
//      break;
//    default:
//      break;
//    }
//    printf("res varval %f oper %d\n",ha_calvar[i].TmpVarVal,ha_calvar[i].Oper);
//  }
//  return ha_calvar[i-1].TmpVarVal;
//}


int ha_newfppow(char *fomulain, ha_cgeset *ha_set,int npow,int ipar,hcge_cof *ha_cof,uvadd ncof, hcge_cof *ha_var,uvadd nvar,uvadd ncofele,hcge_sumcof *sum_cof,int totalsum,hcge_calvars *ha_calvar,int *ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim) {
  int i,i1,i2,i3,i4,i0,ibar=0,index,j,j1,i5,p1;
  char *p=NULL;//,*p1=NULL,*p2=NULL,*p3=NULL,*p4=NULL;
  char fpart1[TABREADLINE],fpart2[TABREADLINE],fpart3[TABREADLINE],var1[TABREADLINE],var2[TABREADLINE],interchar[TABREADLINE],interchar1[TABREADLINE];
  //ha_cgetype var1val=0,var2val=0;

  for (i=1; i<npow+1; i++) {
    p=strchr(fomulain,'^');
    index=p-fomulain;
    strncpy(fpart1, fomulain, index);
    fpart1[index] = '\0';
    strcpy(fpart3,fomulain+index+1);

//    i0=-1;
//    p=strrchr(fpart1,'^');
//    if (p!=NULL) {
//      i0=p-fpart1;
//    }
//    i1=-1;
//    p=strrchr(fpart1,'*');
//    if (p!=NULL) {
//      i1=p-fpart1;
//    }
//    i2=-1;
//    p=strrchr(fpart1,'/');
//    if (p!=NULL) {
//      i2=p-fpart1;
//    }
//    i3=-1;
//    p=strrchr(fpart1,'+');
//    if (p!=NULL) {
//      i3=p-fpart1;
//    }
//    i4=-1;
//    p=strrchr(fpart1,'-');
//    if (p!=NULL) {
//      i4=p-fpart1;
//    }
    i1=-1;
    p=ha_revstrpbrk(fpart1,"^*/+-=<>");
    if (p!=NULL) {
      i1=p-fpart1;
    }
    i5=-1;
    j=strlen(fpart1)-1;
    while(j>-1){
      if(fpart1[j]==','){
        p1=-1;
        for(j1=j-1;j1>-1;j1--){
          if(fpart1[j1]=='{')p1=j1;
          if(fpart1[j1]=='}')break;
        }
        if(p1==-1)break;
        else j=p1;
      }
      j--;
    }
    if(j>0)i5=j;
    if (i1==i5==-1) {//if (i0==i1==i2==i3==i4==i5==-1) {
      index=0;
    } else {
      index=i1;
//      if (index<i2) {
//        index=i2;
//      }
//      if (index<i3) {
//        index=i3;
//      }
//      if (index<i4) {
//        index=i4;
//      }
//      if (index<i0) {
//        index=i0;
//      }
      if (index<i5) {
        index=i5;
      }
    }
    strcpy(var1, fpart1+index+1);
    strncpy(fpart1,fomulain,index+1);
    fpart1[index+1] = '\0';

    ibar=0;
    while (fpart3[ibar] != '\0') {
      ibar++;
    }
//    i0=-1;
//    p=strstr(fpart3,"^");
//    if (p!=NULL) {
//      i0=p-fpart3;
//    }
//    i1=-1;
//    p=strstr(fpart3,"*");
//    if (p!=NULL) {
//      i1=p-fpart3;
//    }
//    i2=-1;
//    p=strstr(fpart3,"/");
//    if (p!=NULL) {
//      i2=p-fpart3;
//    }
//    i3=-1;
//    p=strstr(fpart3,"+");
//    if (p!=NULL) {
//      i3=p-fpart3;
//    }
//    i4=-1;
//    p=strstr(fpart3,"-");
//    if (p!=NULL) {
//      i4=p-fpart3;
//    }
//    if (i0==i1==i2==i3==i4==-1) {
//      index=ibar;
//    } else {
//      if (i1==-1) {
//        index=ibar;
//      } else {
//        index=i1;
//      }
//      if (index>i2&&i2!=-1) {
//        index=i2;
//      }
//      if (index>i3&&i3!=-1) {
//        index=i3;
//      }
//      if (index>i4&&i4!=-1) {
//        index=i4;
//      }
//      if (index>i0&&i0!=-1) {
//        index=i0;
//      }
//    }
    p=strpbrk(fpart3,"^*/+-=<>");
    if(p==NULL)index=ibar;
    else index=p-fpart3;
    strncpy(var2, fpart3, index);
    var2[index] = '\0';
    strcpy(fpart2,fpart3+index);

    hnew_varrepl(var1,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,*ha_calvarsize,arSet,fdim,1);
    hnew_varrepl(var2,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,*ha_calvarsize,arSet,fdim,2);

    sprintf(interchar1, "%d", ipar);
    interchar[0]='\0';
    if (i<10) {
      strcat(interchar,"ha_cgepow");
      strcat(interchar,interchar1);
      strcat(interchar,"000");
    }
    if (9<i&&i<100) {
      strcat(interchar,"ha_cgepow");
      strcat(interchar,interchar1);
      strcat(interchar,"00");
    }
    if (99<i&&i<1000) {
      strcat(interchar,"ha_cgepow");
      strcat(interchar,interchar1);
      strcat(interchar,"0");
    }
    if (999<i&&i<10000) {
      strcat(interchar,"ha_cgepow");
      strcat(interchar,interchar1);
    }
    sprintf(interchar1, "%d", i);
    strcat(interchar,interchar1);
    //printf("interchar %s",interchar);
    strcpy(ha_calvar[*ha_calvarsize].TmpVarName,interchar);
    ha_calvar[*ha_calvarsize].Oper=5;
    //ha_calvar[*ha_calvarsize].TmpVarType=1;
    //ha_calvar[*ha_calvarsize].Var1BegAdd=*ha_calvarsize-1;
    *ha_calvarsize=*ha_calvarsize+1;
    strcat(fpart1, interchar);
    strcat(fpart1, fpart2);
    strcpy(fomulain,fpart1);
  }
  return 1;
}
int ha_newfpmuldiv(char *fomulain, ha_cgeset *ha_set,int nmul,int ipar,hcge_cof *ha_cof,uvadd ncof, hcge_cof *ha_var,uvadd nvar,uvadd ncofele,hcge_sumcof *sum_cof,int totalsum,hcge_calvars *ha_calvar,int *ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim) {
  int i,i1,i2,i3,i4,i5,ibar=0,index,j,j1,p1;
  char *p=NULL;//,*p1=NULL,*p2=NULL,*p3=NULL,*p4=NULL;
  char fpart1[TABREADLINE],fpart2[TABREADLINE],fpart3[TABREADLINE],var1[TABREADLINE],var2[TABREADLINE],interchar[TABREADLINE],interchar1[TABREADLINE];
  //ha_cgetype var1val=0,var2val=0;
  //printf("for %s\n",fomulain);

  for (i=1; i<nmul+1; i++) {
    index=0;
    p=strpbrk(fomulain,"*/");
    if (*p=='/') {
      ha_calvar[*ha_calvarsize].Oper=2;
    } else {
      ha_calvar[*ha_calvarsize].Oper=1;
    }
    index=p-fomulain;

    strncpy(fpart1, fomulain, index);
    fpart1[index] = '\0';
    strcpy(fpart3,fomulain+index+1);
    //printf("fpart1 %s\n",fpart1);
//    i1=-1;
//    p=strrchr(fpart1,'*');
//    if (p!=NULL) {
//      i1=p-fpart1;
//    }
//    i2=-1;
//    p=strrchr(fpart1,'/');
//    if (p!=NULL) {
//      i2=p-fpart1;
//    }
//    i3=-1;
//    p=strrchr(fpart1,'+');
//    if (p!=NULL) {
//      i3=p-fpart1;
//    }
//    i4=-1;
//    p=strrchr(fpart1,'-');
//    if (p!=NULL) {
//      i4=p-fpart1;
//    }
    i1=-1;
    p=ha_revstrpbrk(fpart1,"*/+-=<>");
    if (p!=NULL) {
      i1=p-fpart1;
    }
    i5=-1;
    j=strlen(fpart1)-1;
    while(j>-1){
      if(fpart1[j]==','){
        p1=-1;
        for(j1=j-1;j1>-1;j1--){
          if(fpart1[j1]=='{')p1=j1;
          if(fpart1[j1]=='}')break;
        }
        if(p1==-1)break;
        else j=p1;
      }
      j--;
    }
    if(j>0)i5=j;
    if (i1==i5==-1) {//if (i1==i2==i3==i4==i5==-1) {
      index=0;
    } else {
      index=i1;
//      if (index<i2) {
//        index=i2;
//      }
//      if (index<i3) {
//        index=i3;
//      }
//      if (index<i4) {
//        index=i4;
//      }
      if (index<i5) {
        index=i5;
      }
    }
    strcpy(var1, fpart1+index+1);
    strncpy(fpart1,fomulain,index+1);
    fpart1[index+1] = '\0';

    ibar=0;
    while (fpart3[ibar] != '\0') {
      ibar++;
    }
//    i1=-1;
//    p=strchr(fpart3,'*');
//    if (p!=NULL) {
//      i1=p-fpart3;
//    }
//    i2=-1;
//    p=strchr(fpart3,'/');
//    if (p!=NULL) {
//      i2=p-fpart3;
//    }
//    i3=-1;
//    p=strchr(fpart3,'+');
//    if (p!=NULL) {
//      i3=p-fpart3;
//    }
//    i4=-1;
//    p=strchr(fpart3,'-');
//    if (p!=NULL) {
//      i4=p-fpart3;
//    }
//    if (i1==i2==i3==i4==-1) {
//      index=ibar;
//    } else {
//      if (i1==-1) {
//        index=ibar;
//      } else {
//        index=i1;
//      }
//      if (index>i2&&i2!=-1) {
//        index=i2;
//      }
//      if (index>i3&&i3!=-1) {
//        index=i3;
//      }
//      if (index>i4&&i4!=-1) {
//        index=i4;
//      }
//    }
    p=strpbrk(fpart3,"*/+-=<>");
    if(p==NULL)index=ibar;
    else index=p-fpart3;
    strncpy(var2, fpart3, index);
    var2[index] = '\0';
    strcpy(fpart2,fpart3+index);

    //printf("var1 %s var2 %s\n",var1,var2);
    hnew_varrepl(var1,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,*ha_calvarsize,arSet,fdim,1);
    hnew_varrepl(var2,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,*ha_calvarsize,arSet,fdim,2);

    sprintf(interchar1, "%d", ipar);
    interchar[0]='\0';
    if (i<10) {
      strcat(interchar,"ha_cgemul");
      strcat(interchar,interchar1);
      strcat(interchar,"000");
    }
    if (9<i&&i<100) {
      strcat(interchar,"ha_cgemul");
      strcat(interchar,interchar1);
      strcat(interchar,"00");
    }
    if (99<i&&i<1000) {
      strcat(interchar,"ha_cgemul");
      strcat(interchar,interchar1);
      strcat(interchar,"0");
    }
    if (999<i&&i<10000) {
      strcat(interchar,"ha_cgemul");
      strcat(interchar,interchar1);
    }
    sprintf(interchar1, "%d", i);
    strcat(interchar,interchar1);
    //printf("interchar %s\n",interchar);
    //strcpy(varmul[i-1].varname,interchar);
    strcpy(ha_calvar[*ha_calvarsize].TmpVarName,interchar);
    //ha_calvar[*ha_calvarsize].TmpVarType=1;
    //ha_calvar[*ha_calvarsize].Var1BegAdd=*ha_calvarsize-1;
    *ha_calvarsize=*ha_calvarsize+1;
    strcat(fpart1, interchar);
    strcat(fpart1, fpart2);
    strcpy(fomulain,fpart1);
  }
  return 1;
}

int ha_newfpplumin(char *fomulain, ha_cgeset *ha_set,int nplu,int ipar,hcge_cof *ha_cof,uvadd ncof, hcge_cof *ha_var,uvadd nvar,uvadd ncofele,hcge_sumcof *sum_cof,int totalsum,hcge_calvars *ha_calvar,int *ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim) {
  int i,i1,i3,i4,i5,ibar=0,index,j,j1,p1;
  char *p=NULL;//,*p1=NULL,*p2=NULL,*p3=NULL,*p4=NULL;
  char fpart1[TABREADLINE],fpart2[TABREADLINE],fpart3[TABREADLINE],var1[TABREADLINE],var2[TABREADLINE],interchar[TABREADLINE],interchar1[TABREADLINE];
  //ha_cgetype var1val=0,var2val=0;
  //printf("for %s\n",fomulain);

  for (i=1; i<nplu+1; i++) {
    index=0;
    p=strpbrk(fomulain,"+-");
    if (*p=='+') {
      ha_calvar[*ha_calvarsize].Oper=3;
    } else {
      ha_calvar[*ha_calvarsize].Oper=4;
    }
    index=p-fomulain;

    strncpy(fpart1, fomulain, index);
    fpart1[index] = '\0';
    strcpy(fpart3,fomulain+index+1);
    //printf("fpart1plu %s oper %d\n",fpart1,ha_calvar[*ha_calvarsize].Oper);
//    i3=-1;
//    p=strrchr(fpart1,'+');
//    if (p!=NULL) {
//      i3=p-fpart1;
//    }
//    i4=-1;
//    p=strrchr(fpart1,'-');
//    if (p!=NULL) {
//      i4=p-fpart1;
//    }
    i1=-1;
    p=ha_revstrpbrk(fpart1,"+-=<>");
    if (p!=NULL) {
      i1=p-fpart1;
    }
    i5=-1;
    j=strlen(fpart1)-1;
    while(j>-1){
      if(fpart1[j]==','){
        p1=-1;
        for(j1=j-1;j1>-1;j1--){
          if(fpart1[j1]=='{')p1=j1;
          if(fpart1[j1]=='}')break;
        }
        if(p1==-1)break;
        else j=p1;
      }
      j--;
    }
    if(j>0)i5=j;
    if (i1==i5==-1) {//if (i3==i4==i5==-1) {
      index=0;
    } else {
      index=i1;
//      index=i3;
//      if (index<i4) {
//        index=i4;
//      }
      if (index<i5) {
        index=i5;
      }
    }
    strcpy(var1, fpart1+index+1);
    //printf("indx %d var1 %s fpart3 %s\n",index,var1,fpart3);
    strncpy(fpart1,fomulain,index+1);
    fpart1[index+1] = '\0';

    ibar=0;
    while (fpart3[ibar] != '\0') {
      ibar++;
    }
//    i3=-1;
//    p=strchr(fpart3,'+');
//    if (p!=NULL) {
//      i3=p-fpart3;
//    }
//    i4=-1;
//    p=strchr(fpart3,'-');
//    if (p!=NULL) {
//      i4=p-fpart3;
//    }
//    if (i3==i4==-1) {
//      index=ibar;
//    } else {
//      if (i3==-1) {
//        index=ibar;
//      } else {
//        index=i3;
//      }
//      if (index>i4&&i4!=-1) {
//        index=i4;
//      }
//    }
    p=strpbrk(fpart3,"+-=<>");
    if(p==NULL)index=ibar;
    else index=p-fpart3;
    //printf("var1as %s index %d i3 %d i4 %d ibar %d fpart3 %s\n",var1,index,i3,i4,ibar,fpart3);
    strncpy(var2, fpart3, index);
    var2[index] = '\0';
    //printf("var1 %s var2 %s\n",var1,var2);
    strcpy(fpart2,fpart3+index);

    if(i==1&&var1[0]=='\0'){
      //printf("var1sss %s var2 %s\n",var1,var2);
        //ha_calvar[*ha_calvarsize].Oper=0;
        ha_calvar[*ha_calvarsize].Var1Type=5;
        ha_calvar[*ha_calvarsize].Var1Val=0;
    }
    else hnew_varrepl(var1,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,*ha_calvarsize,arSet,fdim,1);
    if(i==1&&var2[0]=='\0'){
        ha_calvar[*ha_calvarsize].Var2Type=5;
        ha_calvar[*ha_calvarsize].Var2Val=0;
    }
    else hnew_varrepl(var2,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,*ha_calvarsize,arSet,fdim,2);

    sprintf(interchar1, "%d", ipar);
    interchar[0]='\0';
    if (i<10) {
      strcat(interchar,"ha_cgeplu");
      strcat(interchar,interchar1);
      strcat(interchar,"000");
    }
    if (9<i&&i<100) {
      strcat(interchar,"ha_cgeplu");
      strcat(interchar,interchar1);
      strcat(interchar,"00");
    }
    if (99<i&&i<1000) {
      strcat(interchar,"ha_cgeplu");
      strcat(interchar,interchar1);
      strcat(interchar,"0");
    }
    if (999<i&&i<10000) {
      strcat(interchar,"ha_cgeplu");
      strcat(interchar,interchar1);
    }
    sprintf(interchar1, "%d", i);
    strcat(interchar,interchar1);
    //printf("interchar %s\n",interchar);
    //strcpy(varmul[i-1].varname,interchar);
    strcpy(ha_calvar[*ha_calvarsize].TmpVarName,interchar);
    //ha_calvar[*ha_calvarsize].TmpVarType=1;
    //ha_calvar[*ha_calvarsize].Var1BegAdd=*ha_calvarsize-1;
    *ha_calvarsize=*ha_calvarsize+1;
    strcat(fpart1, interchar);
    strcat(fpart1, fpart2);
    strcpy(fomulain,fpart1);
  }
  return 1;
}

int ha_newfpplumin01(char *fomulain, ha_cgeset *ha_set,int nplu,int ipar,hcge_cof *ha_cof,uvadd ncof, hcge_cof *ha_var,uvadd nvar,uvadd ncofele,hcge_sumcof *sum_cof,int totalsum,hcge_calvars *ha_calvar,int *ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim) {
  char *p=NULL,*p1,p2[TABREADLINE],sign='\0',tail[TABREADLINE];
  int i;//,varindex;
  p1=fomulain;
  for (i=0; i<nplu+1; i++) {
    if (i==nplu) {
      strcpy(p2,p1);
    } else {
      p=strpbrk(p1,"+-");
      strncpy(p2,p1,p-p1);
      p2[p-p1]='\0';
      p1=&p[1];
    }
    if(i==0) {
      if(p2[0]=='\0') {
        ha_calvar[*ha_calvarsize].Oper=0;
        ha_calvar[*ha_calvarsize].Var1Type=5;
        ha_calvar[*ha_calvarsize].Var1Val=0;
        ha_calvar[*ha_calvarsize].TmpVarName[0]='\0';
      } else {
        ha_calvar[*ha_calvarsize].Oper=0;
        ha_calvar[*ha_calvarsize].TmpVarName[0]='\0';
        hnew_varrepl(p2,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,*ha_calvarsize,arSet,fdim,1);
      }
    } else {
      if(sign=='+') {
        ha_calvar[*ha_calvarsize].Oper=3;
      } else {
        ha_calvar[*ha_calvarsize].Oper=4;
      }
      ha_calvar[*ha_calvarsize].TmpVarName[0]='\0';
      ha_calvar[*ha_calvarsize].Var1Type=4;
      ha_calvar[*ha_calvarsize].Var1BegAdd=*ha_calvarsize-1;
      //printf("p2a %s\n",p2);
      hnew_varrepl(p2,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,*ha_calvarsize,arSet,fdim,1);
      //printf("p2b %s\n",p2);
    }
    *ha_calvarsize=*ha_calvarsize+1;
    sign=p[0];
  }
  return 1;
}

int ha_newfpif(char *fomulain, ha_cgeset *ha_set,int nif,int ipar,hcge_cof *ha_cof,uvadd ncof, hcge_cof *ha_var,uvadd nvar,uvadd ncofele,hcge_sumcof *sum_cof,int totalsum,hcge_calvars *ha_calvar,int *ha_calvarsize,ha_cgesetindx *arSet,uvdim fdim) {
  char *p=NULL,*p1,*p3,var1[NAMESIZE],var2[NAMESIZE],var3[NAMESIZE];
  int i,j1,j2,j3,l;//,varindex;
  p1=fomulain;
  //printf("v1 %s v2 %s v3 %s\n",fomulain,var2,var3);
  p=strpbrk(p1,"=<>");
  strncpy(var1,p1,p-p1);
  var1[p-p1]='\0';
  //if(*p=='<')printf("p %d\n",*p);
  p3=p+1;
  if(*p3=='='){
    if(*p=='<')ha_calvar[*ha_calvarsize].Oper=75;
    if(*p=='>')ha_calvar[*ha_calvarsize].Oper=76;
    p++;
    p++;
  }else if(*p3=='>'){
          ha_calvar[*ha_calvarsize].Oper=74;
          p++;
          p++;
          }else{
          if(*p=='=')ha_calvar[*ha_calvarsize].Oper=71;
          if(*p=='>')ha_calvar[*ha_calvarsize].Oper=72;
          if(*p=='<')ha_calvar[*ha_calvarsize].Oper=73;
          p++;
        }
//   if (*(p=='='){
//     //ha_calvar[*ha_calvarsize].Oper=71;
//     p3=p-1;
//     if(*p3=='<')ha_calvar[*ha_calvarsize].Oper=75;
//     else if(*p3=='>')ha_calvar[*ha_calvarsize].Oper=76;
//          else ha_calvar[*ha_calvarsize].Oper=71;
//   }
//   if (*p=='>'){
//     p3=p-1;
//     if(*p3=='<')ha_calvar[*ha_calvarsize].Oper=74;
//     else ha_calvar[*ha_calvarsize].Oper=72;
//   }
//   if(*p=='<')ha_calvar[*ha_calvarsize].Oper=73;
//   p++;
  l=strlen(p);
  j1=-1;j2=-1;j3=-1;
  for(i=0;i<l;i++){
    if(*(p+i)==',')j1=i;
    if(*(p+i)=='}')j2=i;
    if(*(p+i)=='{')j3=i;
    if(j3>-1)if(j1>-1&&j1>j2)break;
    else if(j1>-1)break;
  }
  strncpy(var2,p,j1);
  var2[j1]='\0';
  strcpy(var3,p+j1+1);
  //printf("form %s v1 %s v2 %s v3 %s oper %d p %s\n",fomulain,var1,var2,var3,ha_calvar[*ha_calvarsize].Oper,p);
  hnew_varrepl(var1,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,*ha_calvarsize,arSet,fdim,1);
  hnew_varrepl(var3,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,*ha_calvarsize,arSet,fdim,2);
  ha_calvar[*ha_calvarsize].Var3Type=ha_calvar[*ha_calvarsize].Var2Type;
  ha_calvar[*ha_calvarsize].Var3BegAdd=ha_calvar[*ha_calvarsize].Var2BegAdd;
  for(i=0;i<fdim;i++){
  ha_calvar[*ha_calvarsize].Var3leadlag[i]=ha_calvar[*ha_calvarsize].Var2leadlag[i];
  ha_calvar[*ha_calvarsize].Var3SupSet[i]=ha_calvar[*ha_calvarsize].Var3SupSet[i];
  ha_calvar[*ha_calvarsize].Var3SSIndx[i]=ha_calvar[*ha_calvarsize].Var2SSIndx[i];
  ha_calvar[*ha_calvarsize].Var3ADims[i]=ha_calvar[*ha_calvarsize].Var2ADims[i];
  }
  ha_calvar[*ha_calvarsize].Var3Val=ha_calvar[*ha_calvarsize].Var2Val;
  hnew_varrepl(var2,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,*ha_calvarsize,arSet,fdim,2);
  //printf("i0 %d ini0 %lf\n",*ha_calvarsize,ha_calvar[*ha_calvarsize].TmpVarVal);
  //*ha_calvarsize=*ha_calvarsize+1;
  return 1;
}

uvadd hnew_calcff(char *fname, char *commsyntax,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele,bool IsIni) {
  FILE * filehandle;
  char line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],linecopy[TABREADLINE],condvar[MAXVARDIM][NAMESIZE];
  char vname[NAMESIZE],sumsyntax[NAMESIZE],argu[NAMESIZE],tempset[NAMESIZE];
  char *readitem=NULL,*p=NULL,*p1=NULL;
  uvadd i,i1,i3,i4,l,l2=0,j=0,nsumele,dcountdim1[4*MAXVARDIM],ncond,nloops,logioper[MAXVARDIM],logi,logiantidim[MAXVARDIM][MAXVARDIM],logisup[MAXVARDIM][MAXVARDIM],logivarindx[MAXVARDIM],logivartype[MAXVARDIM],index;//m,
  uvdim fdim,dcount,neqsign=0,sup,varsupsetid[MAXVARDIM];
  int ha_calvarsize=0,totalsum,sumcount=1,npow,nmul,ndiv,nplu,nmin,npar,sumindx,b=0;
  uvadd varantidim[MAXVARDIM],varsubset[MAXVARDIM],vararset[MAXVARDIM];
  //bool check;
  ha_cgetype zerodivide=0,cond[MAXVARDIM],eval;
  bool IsFomIni=false,IsDefFomIni=false;
  ha_cgesetindx *arSet1=NULL;
  hcge_calvars *ha_calvar1= NULL;
  //printf("line %s cofvar %f\n",line,ha_cofvar[23].varval);

  filehandle = fopen(fname,"r");
  while (ha_cgertabl1(commsyntax,filehandle,line,ha_cofvar,ha_cof,ncof,&zerodivide,TABREADLINE)) {
    if (strstr(line,"(default=initial)")!=NULL) IsDefFomIni=true;
    if (strstr(line,"(default")==NULL) {
      IsFomIni=IsDefFomIni;
      if(strstr(line, "(initial)")!=NULL) {
        ha_cgefrstr1(line, "(initial)", "");
        IsFomIni=true;
      }
      if(strstr(line, "(always)")!=NULL) {
        ha_cgefrstr1(line, "(always)", "");
        IsFomIni=false;
      }
      //printf("line %s inivar %lf\n",line,ha_cofvar[2678594].varval);
      if(!(IsFomIni&&!IsIni)) {
        ncond=0;
        for (i=0; i<MAXVARDIM; i++)logioper[i]=0;
        ha_cgefrstr1(line, commsyntax, "");
        while (ha_cgefrstr(line," ", ""));
        while (ha_cgefrchr(line, '[', '('));
        while (ha_cgefrchr(line, ']', ')'));
        while (ha_cgefrchr(line, '{', '('));
        while (ha_cgefrchr(line, '}', ')'));
        strcpy(linecopy,line);
        totalsum=ha_cgenfind(line, "sum(");
        neqsign=ha_cgenchf(line, '=');
        readitem = strtok(line,"=");
        for(i=1; i<neqsign; i++)readitem = strtok(NULL,"=");
        readitem = strtok(NULL,";");
        npow=ha_cgenchf(readitem, '^');
        nmul=ha_cgenchf(readitem, '*');
        ndiv=ha_cgenchf(readitem, '/');
        nmul=nmul+ndiv;
        nplu=ha_cgenchf(readitem, '+');
        nmin=ha_cgenchf(readitem, '-');
        nplu=nplu+nmin;

        strcpy(line,linecopy);
        //printf("line %s cofvar %f\n",line,ha_cofvar[23].varval);
        //printf("line %s\n",line);
        //readitem = strtok(line,"=");
        readitem =ha_cgeeqfind(line,'=');//readitem =strrchr(line,'=');
        line[readitem-line]='\0';
        //printf("line1 %s\n",line);
        fdim=ha_cgenchf(line, '(');
        if (fdim==1) {
          fdim=fdim+1;
        }
        ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (fdim+1,sizeof(ha_cgesetindx));

        //strcpy(line,linecopy);
        //readitem = strtok(line,"=");
        //printf("fdim %d line %s\n",fdim,line);
        nloops=1;
        if (fdim==0) {
          readitem = line;
        } else {
          for (i=0; i<fdim-1; i++) {
            if(i==0) {
              readitem = strtok(line,",");
            } else {
              readitem = strtok(NULL,",");
            }
            readitem = strtok(NULL,",");
            //printf("read0 %s\n",readitem);
            strcpy(arSet[i].arIndx,readitem);
            readitem = strtok(NULL,")");
            //printf("read1 %s\n",readitem);
            if(strchr(readitem,':')==NULL)strcpy(tempset,readitem);
            else {
              p = strchr(readitem,':');
              strncpy(tempset,readitem,p-readitem);
              tempset[p-readitem]='\0';
              //printf("ar %s p %s pm %ld\n",tempset,p,p-readitem);
              if(strchr(readitem,'(')!=NULL) {
                fdim--;
                p++;
                strcpy(condvar[i],p);
                strcat(condvar[i],")");
                //printf("convar %s\n",condvar[ncond]);
                //printf("read %s\n",readitem);
                readitem = strtok(NULL,")");
                //printf("read %s\n",readitem);
                if(strstr(readitem,"=")){
                  logioper[i]=1;
                  if(strstr(readitem,"<="))logioper[i]=5;
                  if(strstr(readitem,">="))logioper[i]=6;
                }else{
                  if(strstr(readitem,">"))logioper[i]=2;
                  if(strstr(readitem,"<"))logioper[i]=3;
                  if(strstr(readitem,"<>"))logioper[i]=4;
                }
                //printf("read %s\n",readitem);
                cond[i]=atof(readitem);
                //printf("i %ld convar %s logi %ld\n",i,condvar[i],logioper[i]);
              } else {
                p++;
                p1=strstr(readitem,"=");
                if(p1!=NULL) {
                  logioper[i]=1;
                  strncpy(condvar[i],p,p1-p);
                  cond[i]=atof(readitem);
                p1=strstr(readitem,"<=");
                if(p1!=NULL) {
                  logioper[i]=5;
                  strncpy(condvar[i],p,p1-p);
                  cond[i]=atof(readitem);
                }
                p1=strstr(readitem,">=");
                if(p1!=NULL) {
                  logioper[i]=6;
                  strncpy(condvar[i],p,p1-p);
                  cond[i]=atof(readitem);
                }
                }else{
                p1=strstr(readitem,">");
                if(p1!=NULL) {
                  logioper[i]=2;
                  strncpy(condvar[i],p,p1-p);
                  cond[i]=atof(readitem);
                }
                p1=strstr(readitem,"<");
                if(p1!=NULL) {
                  logioper[i]=3;
                  strncpy(condvar[i],p,p1-p);
                  cond[i]=atof(readitem);
                }
                p1=strstr(readitem,"<>");
                if(p1!=NULL) {
                  logioper[i]=4;
                  strncpy(condvar[i],p,p1-p);
                  cond[i]=atof(readitem);
                }
                }
                //printf("convar %s cond %f logi %ld\n",condvar[i],cond[i],logioper[i]);
              }
              ncond++;
            }
            for (i4=0; i4<nset; i4++) if(strcmp(tempset,ha_set[i4].setname)==0) {
                //printf("set %s size %d\n",ha_set[i4].setname,ha_set[i4].size);
                arSet[i].setid=i4;
                //arSet[i].SetSize=ha_set[i4].size;
                //arSet[i].SetBegAdd=ha_set[i4].begadd;
                //arSet[i].subsetid=ha_set[i4].subsetid;
                //if(ha_set[i4].subsetid==1) {
                //  arSet[i].SuperSetSize=ha_set[i4].supersetsize;
                //} else {
                //  arSet[i].SuperSetSize=ha_set[i4].size;
                //}
                break;
              }
            //printf("set %s n1 %ld ",ha_set[arSet[i].setid].setname,nloops);
            nloops=nloops*ha_set[arSet[i].setid].size;
            //printf("n2 %ld\n",nloops);
          }
          readitem = strtok(NULL,"=");
          //printf("read %s\n",readitem);
          dcountdim1[fdim-2]=1;
          for (i=fdim-3; i>-1; i--) {
            dcountdim1[i]=ha_set[arSet[i+1].setid].size*dcountdim1[i+1];
          }
        }
        //if(ncond>0)printf("ncond %d fdim %d\n",ncond,fdim);
        //if(ncond>0)for (i=0; i<fdim-1; i++)printf("convar %s cond %f logi %d\n",condvar[i],cond[i],logioper[i]);
        //printf("read23 %s\n",readitem);
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        //readitem = strtok(line,"=");
        //for(i=1; i<neqsign; i++)readitem = strtok(NULL,"=");
        //readitem = strtok(line,"=");
        //readitem = strtok(NULL,";");
        readitem =ha_cgeeqfind(line,'=');
        readitem++;
        readitem = strtok(readitem,";");
        //printf("line1 %s\n",readitem);
        while (ha_cgerecovar(readitem)==1);
        //printf("line1 %s\n",readitem);
        hnew_intrpl(readitem);
        //printf("line2 %s\n",readitem);
        npar=ha_cgenchf(readitem, '(');
        strcpy(sumsyntax,"sum(");
        totalsum=hcge_nsum(readitem,sumsyntax);
        hcge_sumcof *sum_cof= (hcge_sumcof *) calloc (totalsum,sizeof(hcge_sumcof));
        sumcount=0;
        strcpy(line1,readitem);
        //printf("line1 %s\n",line1);
        strcpy(line2,line1);
        readitem=line2;
        while (hcge_dsum(readitem,sumsyntax,sum_cof,arSet,ha_set,nset,fdim,sumcount)==1) {
          sumcount++;
        }
        totalsum=sumcount;
        i3=0;
        for (i=0; i<totalsum; i++) {
          i1=1;
          for(j=0; j<sum_cof[i].size; j++) {
            i1=i1*ha_set[sum_cof[i].setid[j]].size;
          }
          sum_cof[i].begadd=i3;
          sum_cof[i].summatsize=i1;
          i3=i3+i1;
        }
        nsumele=i3;
        for (i=0; i<totalsum; i++) {
          i1=1;
          sum_cof[i].antidims[sum_cof[i].size-1]=1;
          for(j=sum_cof[i].size-2; j>-1; j--) {
            sum_cof[i].antidims[j]=sum_cof[i].antidims[j+1]*ha_set[sum_cof[i].setid[j+1]].size;
          }
        }
        hcge_calvars *ha_calvar= (hcge_calvars *) calloc (npow+nmul+nplu+2*(npar+2),sizeof(hcge_calvars));
        ha_cgesumele *ha_sumele= (ha_cgesumele *) calloc (nsumele,sizeof(ha_cgesumele));
        //printf("cs0 %d\n",npow+nmul+nplu+2*npar+2);
        //printf("read2 %s\n",readitem);
        //hcge_rsumele(sum_cof,totalsum,ha_set,nset,ha_setele,ha_sumele);
        //printf("read1 %s\n",readitem);
        sumcount=0;
        strcpy(line2,line1);
        readitem=line2;
        sumindx=0;
        while (hnew_calsum(readitem,sumsyntax,ha_set,nset,ha_setele,ha_cofvar,ncofvar,ncofele,ha_cof,ncof,ha_var,nvar,sum_cof,totalsum,ha_sumele,nsumele,ha_calvar,arSet,fdim,&sumindx,sumcount,zerodivide)==1) {
          sumcount++;
          //if(strcmp(vname,"v1ptx1(c,i79)")==0)
            //printf("vname %s read %s nsumele %d\n",vname,readitem,nsumele);
        }
        //printf("read %s\n",readitem);
        //printf("line %s\n",linecopy);
        //if(strcmp(vname,"v1ptx1(c,i79)")==0)for (i=0;i<nsumele;i++) printf("varval %f\n",ha_sumele[i].varval);
        strcpy(line1,readitem);
        uvadd index=ncof-1, begadd=0;//,simpl=0;
        bool check10=true;
        //if (strpbrk(readitem,"*+-^/)")==NULL) {
        //simpl=1;
        //}
        uvadd varsize=0;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            ha_cof[index].suplval=true;
            begadd=ha_cof[index].begadd;
            varsize=ha_cof[index].size;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            check10=false;
            break;
          }
        } while (index--);
        if (check10) {
          index=nvar-1;
          do {
            if (strcmp(ha_var[index].cofname,p)==0) {
              begadd=ncofele+ha_var[index].begadd;
              varsize=ha_var[index].size;
              ha_var[index].suplval=true;
              if(ha_var[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
              break;
            }
          } while (index--);
        }
        //printf("varindx %ld\n",index);
        //if(strcmp("rgdpcons",p)==0){
        //printf("line %s\n",linecopy);
        //for (i=0;i<totalsum;i++) printf("name %s\n",sum_cof[i].sumname);
        //for (i=0;i<totalsum;i++)for(j=0;j<sum_cof[i].summatsize;j++) printf("varname %s varval %f\n",sum_cof[i].sumname,ha_sumele[sum_cof[i].begadd+j].varval);
        //for (i=0;i<totalsum;i++){printf("name %s ",sum_cof[i].sumname);for(j=0;j<sum_cof[i].size;j++)printf("set %s ",ha_set[sum_cof[i].setid[j]].setname);printf("\n");}
        //}
        //uvadd *varantidim= (uvadd *) calloc (varsize,sizeof(uvadd));
        //uvadd *varsubset= (uvadd *) calloc (varsize,sizeof(uvadd));
        //uvadd *vararset= (uvadd *) calloc (varsize,sizeof(uvadd));
        //printf("argu OK!!!! %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){varantidim[l]=0;varsubset[l]=0;varsupsetid[l]=0;}
        if (check10) {
            for (dcount=0; dcount<ha_var[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  varantidim[dcount]=ha_var[index].antidims[dcount];
                  vararset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_var[index].setid[dcount]].size){
                    //varsubset[dcount]=0;
                  //} else {
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                    //printf("arset %s varsupsetid %d\n",ha_set[arSet[l].setid].setname,varsupsetid[dcount]);
                  }
                  break;
                }
            }
        } else {
            //printf("cof %s argu %s\n",ha_cof[index].cofname,argu);
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  varantidim[dcount]=ha_cof[index].antidims[dcount];
                  vararset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    //varsubset[dcount]=0;
                  //} else {
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                    //printf("arset %s sup0 %d sup1 %s sup2 %d\n",ha_set[arSet[l].setid].setname,ha_set[arSet[l].setid].subsetid[0],ha_set[ha_set[arSet[l].setid].subsetid[1]].setname,ha_set[arSet[l].setid].subsetid[2]);
                    //printf("arset %s supset %d varsupsetid %d cof set size %d\n",ha_set[arSet[l].setid].setname,ha_set[arSet[l].setid].subsetid[sup],varsupsetid[dcount],ha_cof[index].setid[dcount]);
                    //printf("sup %d\n",sup);
                  }
                  break;
                }
            }
        }
        //if (strcmp(vname,"shrem") == 0)printf("vname %s read %s\n line1 %s fdim %d\n\n",vname,readitem,line1,fdim);
        ha_newfparse(line1,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,&ha_calvarsize,arSet,fdim-1);
        //if (strcmp(vname,"shrem") == 0)for(l=0; l<ha_calvarsize; l++) for(l2=0; l2<10; l2++) {
        //printf("name %s oper %d tmptype1 %d begadd1 %d ad1 %d tmptype2 %d begadd2 %d ad2 %d\n",ha_calvar[l].TmpVarName,ha_calvar[l].Oper,ha_calvar[l].Var1Type,ha_calvar[l].Var1BegAdd,ha_calvar[l].Var1ADims[l2],ha_calvar[l].Var2Type,ha_calvar[l].Var2BegAdd,ha_calvar[l].Var2ADims[l2]);
        //}
        if(ncond>0) {
          for(i=0; i<MAXVARDIM; i++)for(j=0; j<MAXVARDIM; j++){
            logiantidim[i][j]=0;
            logisup[i][j]=0;
          }
          //printf("index %ld cof %s size %ld begadd %ld\n",index,ha_cof[index].cofname,ha_cof[index].matsize,ha_cof[index].begadd);
          for(i1=0; i1<fdim-1; i1++) if(logioper[i1]>0){
            index=ncof-1;
            p=strtok(condvar[i1],"(");
            b=0;
            do {
              if (strcmp(ha_cof[index].cofname,condvar[i1])==0) {
                if(!ha_cof[index].suplval)printf("Warning!!!! coefficient %s has not been supplied with values!\n",ha_cof[index].cofname);
                logivarindx[i1]=index;
                logivartype[i1]=0;
                //printf("var2 %s\n",var2);
                b++;
                if(p!=NULL)p=strtok(NULL,")");
                strcpy(argu,p);
                strcat(argu,",");
                for(i=0; i<ha_cof[index].size; i++) {
                  if(i==0)p=strtok(argu,",");
                  else p=strtok(NULL,",");
                  for(j=0; j<fdim; j++) {
                    if(strcmp(arSet[j].arIndx,p)==0) {
                      logiantidim[i1][j]=ha_cof[index].antidims[i];
                      for(sup=1; sup<MAXSUPSET; sup++)if(ha_set[arSet[j].setid].subsetid[sup]==ha_cof[index].setid[i]){logisup[i1][j]=sup;break;}
                      break;
                    }
                  }
                }
                break;
              }
            } while (index--);
            if(b==0) {
              index=nvar-1;
              p=strtok(condvar[i1],"(");
              do {
                if (strcmp(ha_var[index].cofname,condvar[i1])==0) {
                if(!ha_var[index].suplval)printf("Warning!!!! coefficient %s has not been supplied with values!\n",ha_var[index].cofname);
                  logivarindx[i1]=index;
                  logivartype[i1]=1;
                  //printf("var2 %s\n",var2);
                  if(p!=NULL)p=strtok(NULL,")");
                  strcpy(argu,p);
                  strcat(argu,",");
                  for(i=0; i<ha_var[index].size; i++) {
                    if(i==0)p=strtok(argu,",");
                    else p=strtok(NULL,",");
                    for(j=0; j<fdim-1; j++) {
                      if(strcmp(arSet[j].arIndx,p)==0) {
                        logiantidim[i1][j]=ha_var[index].antidims[i];
                        for(sup=1; sup<MAXSUPSET; sup++)if(ha_set[arSet[j].setid].subsetid[sup]==ha_var[index].setid[i]){logisup[i1][j]=sup;break;}
                        break;
                      }
                    }
                  }
                  break;
                }
              } while (index--);
            }
          }
          //printf("OK!!!!!\n");
        #pragma omp parallel private(l,l2,i4,dcount,i3,i1,arSet1,logi,index,eval,ha_calvar1) shared(ha_cofvar,arSet)
        {
        if(omp_get_thread_num()!=0){
          arSet1=realloc(arSet1,(fdim+1)*sizeof(ha_cgesetindx));
          memcpy (arSet1,arSet,(fdim+1)*sizeof(ha_cgesetindx));
          ha_calvar1=realloc(ha_calvar1,(npow+nmul+nplu+2*(npar+2))*sizeof(hcge_calvars));
          memcpy (ha_calvar1,ha_calvar,(npow+nmul+nplu+2*(npar+2))*sizeof(hcge_calvars));
        }else{
          ha_calvar1=ha_calvar;
          arSet1=arSet;
        }
        #pragma omp for
          for (l=0; l<nloops; l++) {
            l2=0;
            i4=l;
            //check=false;
            for (dcount=0; dcount<fdim-1; dcount++) {
              i3=(uvadd) i4/dcountdim1[dcount];
              arSet1[dcount].indx=i3;
              i4=i4-i3*dcountdim1[dcount];
              if(varsize<fdim-1||varsize==fdim-1) {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[arSet1[dcount].setid].begadd+i3].setsh[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                    break;
                  }
                }
              } else {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[arSet1[dcount].setid].begadd+i3].setsh[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                  }
                }
              }
            }
            logi=0;
            index=0;
            for(i1=0; i1<fdim-1; i1++) {
              if(logioper[i1]>0){
              for(i=0; i<fdim-1; i++){
                //index+=arSet1[i].indx*logiantidim[i1][i];
                index+=ha_setele[ha_set[arSet1[i].setid].begadd+arSet1[i].indx].setsh[logisup[i1][i]]*logiantidim[i1][i];
                //if(l==0)printf("arset %d logianti %ld\n",arSet1[i].indx,logiantidim[i1][i]);
              }
              if(logivartype[i1]==0)eval=ha_cofvar[ha_cof[logivarindx[i1]].begadd+index].varval;
              else eval=ha_cofvar[ncofele+ha_var[logivarindx[i1]].begadd+index].varval;
              //if(l<10)printf("i1 %ld convar %s type %ld logoper %ld index %ld begadd %ld size %ld eval %lf cond %lf %ld\n",i1,ha_cof[logivarindx[i1]].cofname,logivartype[i1],logioper[i1],index,ha_cof[logivarindx[i1]].begadd,ha_cof[logivarindx[i1]].matsize,eval,cond[i1],ha_cof[logivarindx[i1]].begadd+index);
              if(logioper[i1]==1)if(eval==cond[i1])logi++;
              if(logioper[i1]==2)if(eval>cond[i1])logi++;
              if(logioper[i1]==3)if(eval<cond[i1])logi++;
              if(logioper[i1]==4)if(eval!=cond[i1])logi++;
              if(logioper[i1]==5)if(eval<=cond[i1])logi++;
              if(logioper[i1]==6)if(eval>=cond[i1])logi++;
              }
            }
            //if(index==125)printf("l2 %ld l %ld logi %ld ncond %ld eval %lf begadd %ld nloop %ld\n",l2,l,logi,ncond,eval,begadd,nloops);
            if(logi==ncond)ha_cofvar[begadd+l2].varval=ha_newfpcal(ha_cofvar,ha_set,ha_setele,ha_sumele,ha_calvar1,ha_calvarsize,arSet1,fdim-1,zerodivide);
            //printf("nloop %d l %d varval %f logi %d\n",nloops,l,ha_cofvar[begadd+l2].varval,ncond);}
          }
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
        } else {
        #pragma omp parallel private(l,l2,i4,dcount,i3,i1,arSet1,ha_calvar1) shared(ha_cofvar,arSet)
        {
        //j8=0;
        if(omp_get_thread_num()!=0){
          arSet1=realloc(arSet1,(fdim+1)*sizeof(ha_cgesetindx));
          memcpy (arSet1,arSet,(fdim+1)*sizeof(ha_cgesetindx));
          ha_calvar1=realloc(ha_calvar1,(npow+nmul+nplu+2*(npar+2))*sizeof(hcge_calvars));
          memcpy (ha_calvar1,ha_calvar,(npow+nmul+nplu+2*(npar+2))*sizeof(hcge_calvars));
        }else{
          ha_calvar1=ha_calvar;
          arSet1=arSet;
        }
        #pragma omp for
          for (l=0; l<nloops; l++) {
            l2=0;
            i4=l;
            //check=false;
            for (dcount=0; dcount<fdim-1; dcount++) {
              i3=(uvadd) i4/dcountdim1[dcount];
              arSet1[dcount].indx=i3;
              //printf("i3 %d l2 %d\n",i3,l2);
              i4=i4-i3*dcountdim1[dcount];
              if(varsize<=fdim-1) {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[arSet1[dcount].setid].begadd+i3].setsh[varsupsetid[i1]]*varantidim[i1];
                      //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                    break;
                  }
                }
              } else {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[arSet1[dcount].setid].begadd+i3].setsh[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                  }
                }
              }
            }
            ha_cofvar[begadd+l2].varval=ha_newfpcal(ha_cofvar,ha_set,ha_setele,ha_sumele,ha_calvar1,ha_calvarsize,arSet1,fdim-1,zerodivide);
            //printf("vname %s\n",vname);
            //if (strcmp(vname,"v1prim") == 0)printf("vname %s val %lf\n",vname,ha_cofvar[begadd+l2].varval);
            //for(i1=0;i1<ncof;i1++)if(ha_cof[i1].begadd==begadd)printf("name %s nloops %d l2 %d varval1 %f\n",ha_cof[i1].cofname,nloops,l2,ha_cofvar[begadd+l2].varval);
          }
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
        }

        if(ha_cof[index].gltype>0){
        #pragma omp parallel private(l) shared(ha_cofvar,ha_cof,ncof,begadd,index)
        {
          if(ha_cof[index].gltype==1){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval<ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 1!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
          if(ha_cof[index].gltype==2){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval<=ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 2!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
          if(ha_cof[index].gltype==3){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval>ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 3!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
          if(ha_cof[index].gltype==4){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval>=ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 4!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
        }
        }
        
        
        
        free(sum_cof);
        free(ha_sumele);
        //free(varantidim);
        //free(varsubset);
        //free(vararset);
        free(arSet);
        free(ha_calvar);
      }
    }
  }
  fclose(filehandle);
  return j;
}

uvadd hnew_update(char *fname,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele) {
  FILE * filehandle;
  char commsyntax[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],linecopy[TABREADLINE];
  char vname[NAMESIZE],sumsyntax[NAMESIZE],argu[NAMESIZE];
  char *readitem=NULL,*p=NULL;
  uvadd i,i1,i3,i4,l,l2=0,j=0,nsumele,dcountdim1[4*MAXVARDIM],nloops;//m,
  uvdim fdim,dcount,sup,varsupsetid[MAXVARDIM];
  int ha_calvarsize=0,totalsum,sumcount=1,npow,nmul,ndiv,nplu,nmin,npar,sumindx;
  bool IsChange=false,IsExplicit=false;
  ha_cgetype zerodivide=0,temp1;
  uvadd varantidim[MAXVARDIM],varsubset[MAXVARDIM],vararset[MAXVARDIM];
  //printf("line %s cofvar %f\n",line,ha_cofvar[2083].varval);
  ha_cgesetindx *arSet1=NULL;
  hcge_calvars *ha_calvar1= NULL;
  strcpy(commsyntax,"update");
  filehandle = fopen(fname,"r");
  while (ha_cgertabl1(commsyntax,filehandle,line,ha_cofvar,ha_cof,ncof,&zerodivide,TABREADLINE)) {
    IsChange=false;
    IsExplicit=false;
    if(strstr(line, "(change)")!=NULL) {
      IsChange=true;
      ha_cgefrstr1(line, "(change)", "");
    }
    if(strstr(line, "(explicit)")!=NULL) {
      IsExplicit=true;
      ha_cgefrstr1(line, "(explicit)", "");
    }
    ha_cgefrstr1(line, commsyntax, "");
    while (ha_cgefrstr(line," ", ""));
    while (ha_cgefrchr(line, '[', '('));
    while (ha_cgefrchr(line, ']', ')'));
    while (ha_cgefrchr(line, '{', '('));
    while (ha_cgefrchr(line, '}', ')'));
    strcpy(linecopy,line);
    //printf("line %s cofvar %f\n",line,ha_cofvar[23].varval);
    //printf("line %s\n",line);
    readitem = strtok(line,"=");
    fdim=ha_cgenchf(readitem, '(');
    if (fdim==1) {
      fdim=fdim+1;
    }
    ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (fdim+1,sizeof(ha_cgesetindx));

    nloops=1;
    if (fdim==0) {
      readitem = strtok(line,"=");
    } else {
      for (i=0; i<fdim-1; i++) {
        if(i==0) {
          readitem = strtok(line,",");
        } else {
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
            //} else {
            //  arSet[i].SuperSetSize=ha_set[i4].size;
            //}
            break;
          }
        nloops=nloops*ha_set[arSet[i].setid].size;
      }
      readitem = strtok(NULL,"=");
      dcountdim1[fdim-2]=1;
      for (i=fdim-3; i>-1; i--) {
        dcountdim1[i]=ha_set[arSet[i+1].setid].size*dcountdim1[i+1];
      }
    }
    strcpy(vname,readitem);
    //if(strcmp(vname,"vkb(r)")==0)printf("vname %s\n",vname);
    strcpy(line,linecopy);
    readitem = strtok(line,"=");
    readitem = strtok(NULL,";");
    if(IsChange==false&&IsExplicit==false) {
      //if(lastupd){
      //while (ha_cgefrstr(readitem,"*", "/100)!(1+"));
      //while (ha_cgefrstr(readitem,"!", "*"));
      //strcpy(line1,vname);
      //strcat(line1,"*(1+");
      //strcat(line1,readitem);
      //strcat(line1,"/100)");
      //readitem=line1;
      //} else{
      while (ha_cgefrstr(readitem,"*", "+"));
      strcpy(line1,vname);
      strcat(line1,"*(1+(");
      strcat(line1,readitem);
      strcat(line1,")/100)");
      readitem=line1;
      //}
    }
    if(IsChange==true) {
      strcpy(line1,vname);
      if(readitem[0]!='+'||readitem[0]!='-') strcat(line1,"+");
      strcat(line1,readitem);
      readitem=line1;
    }
    //if(IsExplicit==true){
    //  strcpy(line1,vname);
    //  if(readitem[0]!='+'||readitem[0]!='-') strcat(line1,"+");
    //  strcat(line1,readitem);
    //  readitem=line1;
    //}
    npow=ha_cgenchf(readitem, '^');
    nmul=ha_cgenchf(readitem, '*');
    ndiv=ha_cgenchf(readitem, '/');
    nmul=nmul+ndiv;
    nplu=ha_cgenchf(readitem, '+');
    nmin=ha_cgenchf(readitem, '-');
    nplu=nplu+nmin;

    //strcpy(line,linecopy);
    //readitem = strtok(line,"=");
    //readitem = strtok(NULL,";");
    while (ha_cgerecovar(readitem)==1);
    //printf("line1 %s\n",readitem);
    hnew_intrpl(readitem);
    //printf("line2 %s\n",readitem);
    npar=ha_cgenchf(readitem, '(');
    strcpy(sumsyntax,"sum(");
    totalsum=hcge_nsum(readitem,sumsyntax);
    hcge_sumcof *sum_cof= (hcge_sumcof *) calloc (totalsum,sizeof(hcge_sumcof));
    sumcount=0;
    strcpy(line1,readitem);
    //printf("line1 %s\n",line1);
    strcpy(line2,line1);
    readitem=line2;
    while (hcge_dsum(readitem,sumsyntax,sum_cof,arSet,ha_set,nset,fdim,sumcount)==1) {
      sumcount++;
    }
    totalsum=sumcount;
    i3=0;
    for (i=0; i<totalsum; i++) {
      i1=1;
      for(j=0; j<sum_cof[i].size; j++) {
        i1=i1*ha_set[sum_cof[i].setid[j]].size;
      }
      sum_cof[i].begadd=i3;
      i3=i3+i1;
    }
    nsumele=i3;
    for (i=0; i<totalsum; i++) {
      i1=1;
      sum_cof[i].antidims[sum_cof[i].size-1]=1;
      for(j=sum_cof[i].size-2; j>-1; j--) {
        sum_cof[i].antidims[j]=sum_cof[i].antidims[j+1]*ha_set[sum_cof[i].setid[j+1]].size;
      }
    }
    hcge_calvars *ha_calvar= (hcge_calvars *) calloc (npow+nmul+nplu+2*npar+2,sizeof(hcge_calvars));
    ha_cgesumele *ha_sumele= (ha_cgesumele *) calloc (nsumele,sizeof(ha_cgesumele));
    //printf("read2 %s\n",readitem);
    //hcge_rsumele(sum_cof,totalsum,ha_set,nset,ha_setele,ha_sumele);
    //printf("read1 %s\n",readitem);
    sumcount=0;
    strcpy(line2,line1);
    readitem=line2;
    sumindx=0;
    while (hnew_calsum(readitem,sumsyntax,ha_set,nset,ha_setele,ha_cofvar,ncofvar,ncofele,ha_cof,ncof,ha_var,nvar,sum_cof,totalsum,ha_sumele,nsumele,ha_calvar,arSet,fdim,&sumindx,sumcount,zerodivide)==1) {
      sumcount++;
    }
    //printf("read0 %s\n",readitem);
    //for (i=0;i<nsumele;i++) printf("varele %s varval%f\n",ha_sumele[i].varname,ha_sumele[i].varval);
    strcpy(line1,readitem);
    uvadd index=ncof-1, begadd=0;//,simpl=0;
    bool check10=true;
    //if (strpbrk(readitem,"*+-^/)")==NULL) {
    //simpl=1;
    //}
        uvadd varsize=0;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            begadd=ha_cof[index].begadd;
            varsize=ha_cof[index].size;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            check10=false;
            break;
          }
        } while (index--);
        if (check10) {
          index=nvar-1;
          do {
            if (strcmp(ha_var[index].cofname,p)==0) {
              begadd=ncofele+ha_var[index].begadd;
              varsize=ha_var[index].size;
              if(ha_var[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
              break;
            }
          } while (index--);
        }
        for (l=0; l<MAXVARDIM; l++){varantidim[l]=0;varsubset[l]=0;varsupsetid[l]=0;}
        //printf("OKK check10 %d s %ld p %s\n",check10,index,p);
        if (check10) {
            for (dcount=0; dcount<ha_var[index].size; dcount++) {
              //printf("dcount %d\n",dcount);
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  varantidim[dcount]=ha_var[index].antidims[dcount];
                  vararset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_var[index].setid[dcount]].size){
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
        } else {
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  varantidim[dcount]=ha_cof[index].antidims[dcount];
                  vararset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
        }
    /*uvadd varsize=0;
    p=strtok(vname,"(");
    do {
      if (strcmp(ha_cof[index].cofname,p)==0) {
        begadd=ha_cof[index].begadd;
        varsize=ha_cof[index].size;
        check10=false;
        break;
      }
    } while (index--);
    if (check10) {
      index=nvar-1;
      do {
        if (strcmp(ha_var[index].cofname,p)==0) {
          begadd=ncofele+ha_var[index].begadd;
          varsize=ha_var[index].size;
          break;
        }
      } while (index--);
    }
    //for (l=0; l<MAXVARDIM; l++)varsupid[l]=0;
    if (check10) {
      for (dcount=0; dcount<ha_var[index].size-1; dcount++) {
        p=strtok(NULL,",");
        for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
            varantidim[dcount]=ha_var[index].antidims[dcount];
            vararset[dcount]=l+1;
            if (strcmp(ha_set[arSet[l].setid].setname,ha_set[ha_var[index].setid[dcount]].setname)==0) {
              varsubset[dcount]=0;
            } else {
              varsubset[dcount]=1;
              for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
            }
            break;
          }
      }
      p=strtok(NULL,")");
      for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
          varantidim[dcount]=ha_var[index].antidims[dcount];
          vararset[dcount]=l+1;
          if (strcmp(ha_set[arSet[l].setid].setname,ha_set[ha_var[index].setid[dcount]].setname)==0) {
            varsubset[dcount]=0;
          } else {
            varsubset[dcount]=1;
            for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
          }
          break;
        }
    }else{
      for (dcount=0; dcount<ha_cof[index].size-1; dcount++) {
        p=strtok(NULL,",");
        for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
            varantidim[dcount]=ha_cof[index].antidims[dcount];
            vararset[dcount]=l+1;
            if (strcmp(ha_set[arSet[l].setid].setname,ha_set[ha_cof[index].setid[dcount]].setname)==0) {
              varsubset[dcount]=0;
            } else {
              varsubset[dcount]=1;
              for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
            }
            break;
          }
      }
      p=strtok(NULL,")");
      for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
          varantidim[dcount]=ha_cof[index].antidims[dcount];
          vararset[dcount]=l+1;
          if (strcmp(ha_set[arSet[l].setid].setname,ha_set[ha_cof[index].setid[dcount]].setname)==0) {
            varsubset[dcount]=0;
          } else {
            varsubset[dcount]=1;
            for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
          }
          break;
        }
    }*/
    //printf("read %s\n line1 %s fdim %d\n\n",readitem,line1,fdim);
    ha_newfparse(line1,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,&ha_calvarsize,arSet,fdim-1);
    //printf("clsize %d\n",ha_calvarsize);
    //for(l=0; l<ha_calvarsize; l++) for(l2=0; l2<10; l2++) {
    //printf("name %s oper %d tmptype1 %d begadd1 %d ad1 %d tmptype2 %d begadd2 %d ad2 %d\n",ha_calvar[l].TmpVarName,ha_calvar[l].Oper,ha_calvar[l].Var1Type,ha_calvar[l].Var1BegAdd,ha_calvar[l].Var1ADims[l2],ha_calvar[l].Var2Type,ha_calvar[l].Var2BegAdd,ha_calvar[l].Var2ADims[l2]);
    //}
    //printf("vname %s\n",vname);
        #pragma omp parallel private(l,l2,i4,dcount,i3,i1,temp1,arSet1,ha_calvar1) shared(ha_cofvar,arSet)
        {
        if(omp_get_thread_num()!=0){
          arSet1=realloc(arSet1,(fdim+1)*sizeof(ha_cgesetindx));
          memcpy (arSet1,arSet,(fdim+1)*sizeof(ha_cgesetindx));
          ha_calvar1=realloc(ha_calvar1,ha_calvarsize*sizeof(hcge_calvars));
          memcpy (ha_calvar1,ha_calvar,ha_calvarsize*sizeof(hcge_calvars));
        }else{
          ha_calvar1=ha_calvar;
          arSet1=arSet;
        }
        #pragma omp for
    for (l=0; l<nloops; l++) {
      l2=0;
      i4=l;
      //check=false;
            for (dcount=0; dcount<fdim-1; dcount++) {
              i3=(uvadd) i4/dcountdim1[dcount];
              arSet1[dcount].indx=i3;
              //printf("i3 %d l2 %d\n",i3,l2);
              i4=i4-i3*dcountdim1[dcount];
              if(varsize<=fdim-1) {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[arSet1[dcount].setid].begadd+i3].setsh[varsupsetid[i1]]*varantidim[i1];
                      //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                    break;
                  }
                }
              } else {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[arSet1[dcount].setid].begadd+i3].setsh[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                  }
                }
              }
            }
      /*for (dcount=0; dcount<fdim-1; dcount++) {
        i3=(uvadd) i4/dcountdim1[dcount];
        arSet[dcount].indx=i3;
        i4=i4-i3*dcountdim1[dcount];
        for(i1=0; i1<varsize; i1++) {
          if(vararset[i1]-1==dcount) {
            if(varsubset[i1]==1) {
              l2=l2+ha_setele[ha_set[arSet[dcount].setid].begadd+i3].setsh[varsupid[i1]]*varantidim[i1];
            } else {
              l2=l2+i3*varantidim[i1];
            }
            break;
          }
        }
      }*/
      //printf("varval %f i %d l2 %d\n",ha_cofvar[begadd+l2].varval,begadd+l2,l2);
      //printf("l2 %d var %lf p %lf kb %lf\n",l2,ha_cofvar[begadd+l2].varval,ha_cofvar[813046+106689+l2].csolpupd,ha_cofvar[813046+106656+l2].csolpupd);
      ha_cofvar[begadd+l2].csolpupd=ha_cofvar[begadd+l2].varval;
      temp1=ha_newfpcal(ha_cofvar,ha_set,ha_setele,ha_sumele,ha_calvar1,ha_calvarsize,arSet1,fdim-1,zerodivide);
      if(temp1-ha_cofvar[begadd+l2].varval>0.000000001||temp1-ha_cofvar[begadd+l2].varval<-0.000000001)ha_cofvar[begadd+l2].varval=temp1;
      //printf("varval %f i %d l2 %d old valf %f\n",ha_cofvar[begadd+l2].varval,begadd+l2,l2,temp1);
      //if(strcmp(vname,"vkb")==0)printf("varval %f i %d l2 %d past %lf\n",ha_cofvar[begadd+l2].varval,begadd+l2,l2,ha_cofvar[begadd+l2].csolpupd);
    }
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
    free(sum_cof);
    free(ha_sumele);
    free(arSet);
    free(ha_calvar);

        if(ha_cof[index].gltype>0){
        #pragma omp parallel private(l) shared(ha_cofvar,ha_cof,ncof,begadd,index)
        {
          if(ha_cof[index].gltype==1){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval<ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 1!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
          if(ha_cof[index].gltype==2){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval<=ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 2!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
          if(ha_cof[index].gltype==3){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval>ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 3!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
          if(ha_cof[index].gltype==4){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval>=ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 4!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
        }
        }
    
  }
  fclose(filehandle);
  return j;
}

uvadd hnew_mupdate(char *fname,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele) {
  FILE * filehandle;
  char commsyntax[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],linecopy[TABREADLINE];
  char vname[NAMESIZE],sumsyntax[NAMESIZE],argu[NAMESIZE];
  char *readitem=NULL,*p=NULL;
  uvadd i,i1,i3,i4,l,l2=0,j=0,nsumele,dcountdim1[4*MAXVARDIM],nloops;//m,
  uvdim fdim,dcount,sup,varsupsetid[MAXVARDIM];
  int ha_calvarsize=0,totalsum,sumcount=1,npow,nmul,ndiv,nplu,nmin,npar,sumindx;
  bool IsChange=false,IsExplicit=false;
  ha_cgetype zerodivide=0,temp1,temp2;
  uvadd varantidim[MAXVARDIM],varsubset[MAXVARDIM],vararset[MAXVARDIM];
  //printf("line %s cofvar %f\n",line,ha_cofvar[23].varval);
  ha_cgesetindx *arSet1=NULL;
  hcge_calvars *ha_calvar1= NULL;
  strcpy(commsyntax,"update");
  filehandle = fopen(fname,"r");
  while (ha_cgertabl1(commsyntax,filehandle,line,ha_cofvar,ha_cof,ncof,&zerodivide,TABREADLINE)) {
    IsChange=false;
    IsExplicit=false;
    if(strstr(line, "(change)")!=NULL) {
      IsChange=true;
      ha_cgefrstr1(line, "(change)", "");
    }
    if(strstr(line, "(explicit)")!=NULL) {
      IsExplicit=true;
      ha_cgefrstr1(line, "(explicit)", "");
    }
    ha_cgefrstr1(line, commsyntax, "");
    while (ha_cgefrstr(line," ", ""));
    while (ha_cgefrchr(line, '[', '('));
    while (ha_cgefrchr(line, ']', ')'));
    while (ha_cgefrchr(line, '{', '('));
    while (ha_cgefrchr(line, '}', ')'));
    strcpy(linecopy,line);
    //printf("line %s cofvar %f\n",line,ha_cofvar[23].varval);
    //printf("line %s\n",line);
    readitem = strtok(line,"=");
    fdim=ha_cgenchf(readitem, '(');
    if (fdim==1) {
      fdim=fdim+1;
    }
    ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (fdim+1,sizeof(ha_cgesetindx));

    nloops=1;
    if (fdim==0) {
      readitem = strtok(line,"=");
    } else {
      for (i=0; i<fdim-1; i++) {
        if(i==0) {
          readitem = strtok(line,",");
        } else {
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
            //} else {
            //  arSet[i].SuperSetSize=ha_set[i4].size;
            //}
            break;
          }
        nloops=nloops*ha_set[arSet[i].setid].size;
      }
      readitem = strtok(NULL,"=");
      dcountdim1[fdim-2]=1;
      for (i=fdim-3; i>-1; i--) {
        dcountdim1[i]=ha_set[arSet[i+1].setid].size*dcountdim1[i+1];
      }
    }
    strcpy(vname,readitem);
    //printf("vname %s\n",vname);
    strcpy(line,linecopy);
    readitem = strtok(line,"=");
    readitem = strtok(NULL,";");
    if(IsChange==false&&IsExplicit==false) {
      //if(lastupd){
      //while (ha_cgefrstr(readitem,"*", "/100)!(1+"));
      //while (ha_cgefrstr(readitem,"!", "*"));
      //strcpy(line1,vname);
      //strcat(line1,"*(1+");
      //strcat(line1,readitem);
      //strcat(line1,"/100)");
      //readitem=line1;
      //} else{
      while (ha_cgefrstr(readitem,"*", "+"));
      strcpy(line1,vname);
      strcat(line1,"*(1+(");
      strcat(line1,readitem);
      strcat(line1,")/100)");
      readitem=line1;
      //}
    }
    if(IsChange==true) {
      strcpy(line1,vname);
      if(readitem[0]!='+'||readitem[0]!='-') strcat(line1,"+");
      strcat(line1,readitem);
      readitem=line1;
    }
    npow=ha_cgenchf(readitem, '^');
    nmul=ha_cgenchf(readitem, '*');
    ndiv=ha_cgenchf(readitem, '/');
    nmul=nmul+ndiv;
    nplu=ha_cgenchf(readitem, '+');
    nmin=ha_cgenchf(readitem, '-');
    nplu=nplu+nmin;

    //strcpy(line,linecopy);
    //readitem = strtok(line,"=");
    //readitem = strtok(NULL,";");
    while (ha_cgerecovar(readitem)==1);
    //printf("line1 %s\n",readitem);
    hnew_intrpl(readitem);
    //printf("line2 %s\n",readitem);
    npar=ha_cgenchf(readitem, '(');
    strcpy(sumsyntax,"sum(");
    totalsum=hcge_nsum(readitem,sumsyntax);
    hcge_sumcof *sum_cof= (hcge_sumcof *) calloc (totalsum,sizeof(hcge_sumcof));
    sumcount=0;
    strcpy(line1,readitem);
    //printf("line1 %s\n",line1);
    strcpy(line2,line1);
    readitem=line2;
    while (hcge_dsum(readitem,sumsyntax,sum_cof,arSet,ha_set,nset,fdim,sumcount)==1) {
      sumcount++;
    }
    totalsum=sumcount;
    i3=0;
    for (i=0; i<totalsum; i++) {
      i1=1;
      for(j=0; j<sum_cof[i].size; j++) {
        i1=i1*ha_set[sum_cof[i].setid[j]].size;
      }
      sum_cof[i].begadd=i3;
      i3=i3+i1;
    }
    nsumele=i3;
    for (i=0; i<totalsum; i++) {
      i1=1;
      sum_cof[i].antidims[sum_cof[i].size-1]=1;
      for(j=sum_cof[i].size-2; j>-1; j--) {
        sum_cof[i].antidims[j]=sum_cof[i].antidims[j+1]*ha_set[sum_cof[i].setid[j+1]].size;
      }
    }
    hcge_calvars *ha_calvar= (hcge_calvars *) calloc (npow+nmul+nplu+2*npar+2,sizeof(hcge_calvars));
    ha_cgesumele *ha_sumele= (ha_cgesumele *) calloc (nsumele,sizeof(ha_cgesumele));
    //printf("read2 %s\n",readitem);
    //hcge_rsumele(sum_cof,totalsum,ha_set,nset,ha_setele,ha_sumele);
    //printf("read1 %s\n",readitem);
    sumcount=0;
    strcpy(line2,line1);
    readitem=line2;
    sumindx=0;
    while (hnew_calsum(readitem,sumsyntax,ha_set,nset,ha_setele,ha_cofvar,ncofvar,ncofele,ha_cof,ncof,ha_var,nvar,sum_cof,totalsum,ha_sumele,nsumele,ha_calvar,arSet,fdim,&sumindx,sumcount,zerodivide)==1) {
      sumcount++;
    }
    //printf("read %s\n",readitem);
    //for (i=0;i<nsumele;i++) printf("varele %s varval%f\n",ha_sumele[i].varname,ha_sumele[i].varval);
    strcpy(line1,readitem);
    uvadd index=ncof-1, begadd=0;//,simpl=0;
    bool check10=true;
    //if (strpbrk(readitem,"*+-^/)")==NULL) {
    //simpl=1;
    //}
        uvadd varsize=0;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            begadd=ha_cof[index].begadd;
            varsize=ha_cof[index].size;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            check10=false;
            break;
          }
        } while (index--);
        if (check10) {
          index=nvar-1;
          do {
            if (strcmp(ha_var[index].cofname,p)==0) {
              begadd=ncofele+ha_var[index].begadd;
              varsize=ha_var[index].size;
              if(ha_var[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
              break;
            }
          } while (index--);
        }
        //if(strcmp("rgdpcons",p)==0){
        //printf("line %s\n",linecopy);
        //for (i=0;i<totalsum;i++) printf("name %s\n",sum_cof[i].sumname);
        //for (i=0;i<totalsum;i++)for(j=0;j<sum_cof[i].summatsize;j++) printf("varname %s varval %f\n",sum_cof[i].sumname,ha_sumele[sum_cof[i].begadd+j].varval);
        //for (i=0;i<totalsum;i++){printf("name %s ",sum_cof[i].sumname);for(j=0;j<sum_cof[i].size;j++)printf("set %s ",ha_set[sum_cof[i].setid[j]].setname);printf("\n");}
        //}
        //uvadd *varantidim= (uvadd *) calloc (varsize,sizeof(uvadd));
        //uvadd *varsubset= (uvadd *) calloc (varsize,sizeof(uvadd));
        //uvadd *vararset= (uvadd *) calloc (varsize,sizeof(uvadd));
        //printf("argu OK!!!! %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){varantidim[l]=0;varsubset[l]=0;varsupsetid[l]=0;}
        if (check10) {
            for (dcount=0; dcount<ha_var[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  varantidim[dcount]=ha_var[index].antidims[dcount];
                  vararset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_var[index].setid[dcount]].size){
                    //varsubset[dcount]=0;
                  //} else {
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                    //printf("arset %s varsupsetid %d\n",ha_set[arSet[l].setid].setname,varsupsetid[dcount]);
                  }
                  break;
                }
            }
        } else {
            //printf("cof %s argu %s\n",ha_cof[index].cofname,argu);
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  varantidim[dcount]=ha_cof[index].antidims[dcount];
                  vararset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    //varsubset[dcount]=0;
                  //} else {
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                    //printf("arset %s sup0 %d sup1 %s sup2 %d\n",ha_set[arSet[l].setid].setname,ha_set[arSet[l].setid].subsetid[0],ha_set[ha_set[arSet[l].setid].subsetid[1]].setname,ha_set[arSet[l].setid].subsetid[2]);
                    //printf("arset %s supset %d varsupsetid %d cof set size %d\n",ha_set[arSet[l].setid].setname,ha_set[arSet[l].setid].subsetid[sup],varsupsetid[dcount],ha_cof[index].setid[dcount]);
                    //printf("sup %d\n",sup);
                  }
                  break;
                }
            }
        }
    /*uvadd varsize=0;
    p=strtok(vname,"(");
    do {
      if (strcmp(ha_cof[index].cofname,p)==0) {
        begadd=ha_cof[index].begadd;
        varsize=ha_cof[index].size;
        check10=false;
        break;
      }
    } while (index--);
    if (check10) {
      index=nvar-1;
      do {
        if (strcmp(ha_var[index].cofname,p)==0) {
          begadd=ncofele+ha_var[index].begadd;
          varsize=ha_var[index].size;
          break;
        }
      } while (index--);
    }
    for (l=0; l<MAXVARDIM; l++)varsupid[l]=0;
    if (check10&&ha_var[index].size>0) {
      for (dcount=0; dcount<ha_var[index].size-1; dcount++) {
        p=strtok(NULL,",");
        for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
            varantidim[dcount]=ha_var[index].antidims[dcount];
            vararset[dcount]=l+1;
            if (strcmp(ha_set[arSet[l].setid].setname,ha_set[ha_var[index].setid[dcount]].setname)==0) {
              varsubset[dcount]=0;
            } else {
              varsubset[dcount]=1;
              for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupid[dcount]=sup;break;}
            }
            break;
          }
      }
      p=strtok(NULL,")");
      for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
          varantidim[dcount]=ha_var[index].antidims[dcount];
          vararset[dcount]=l+1;
          if (strcmp(ha_set[arSet[l].setid].setname,ha_set[ha_var[index].setid[dcount]].setname)==0) {
            varsubset[dcount]=0;
          } else {
            varsubset[dcount]=1;
            for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupid[dcount]=sup;break;}
          }
          break;
        }
    }
    if (!check10&&ha_cof[index].size>0) {
      for (dcount=0; dcount<ha_cof[index].size-1; dcount++) {
        p=strtok(NULL,",");
        for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
            varantidim[dcount]=ha_cof[index].antidims[dcount];
            vararset[dcount]=l+1;
            if (strcmp(ha_set[arSet[l].setid].setname,ha_set[ha_cof[index].setid[dcount]].setname)==0) {
              varsubset[dcount]=0;
            } else {
              varsubset[dcount]=1;
              for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupid[dcount]=sup;break;}
            }
            break;
          }
      }
      p=strtok(NULL,")");
      for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
          varantidim[dcount]=ha_cof[index].antidims[dcount];
          vararset[dcount]=l+1;
          if (strcmp(ha_set[arSet[l].setid].setname,ha_set[ha_cof[index].setid[dcount]].setname)==0) {
            varsubset[dcount]=0;
          } else {
            varsubset[dcount]=1;
            for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupid[dcount]=sup;break;}
          }
          break;
        }
    }*/
    //printf("read %s\n line1 %s fdim %d\n\n",readitem,line1,fdim);
    ha_newfparse(line1,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,&ha_calvarsize,arSet,fdim-1);
    //for(l=0; l<ha_calvarsize; l++) for(l2=0; l2<10; l2++) {
    //printf("name %s oper %d tmptype1 %d begadd1 %d ad1 %d tmptype2 %d begadd2 %d ad2 %d\n",ha_calvar[l].TmpVarName,ha_calvar[l].Oper,ha_calvar[l].Var1Type,ha_calvar[l].Var1BegAdd,ha_calvar[l].Var1ADims[l2],ha_calvar[l].Var2Type,ha_calvar[l].Var2BegAdd,ha_calvar[l].Var2ADims[l2]);
    //}
        #pragma omp parallel private(l,l2,i4,dcount,i3,i1,temp1,temp2,arSet1,ha_calvar1) shared(ha_cofvar,arSet)
        {
        if(omp_get_thread_num()!=0){
          arSet1=realloc(arSet1,(fdim+1)*sizeof(ha_cgesetindx));
          memcpy (arSet1,arSet,(fdim+1)*sizeof(ha_cgesetindx));
          ha_calvar1=realloc(ha_calvar1,ha_calvarsize*sizeof(hcge_calvars));
          memcpy (ha_calvar1,ha_calvar,ha_calvarsize*sizeof(hcge_calvars));
        }else{
          ha_calvar1=ha_calvar;
          arSet1=arSet;
        }
        #pragma omp for
    for (l=0; l<nloops; l++) {
      l2=0;
      i4=l;
      //check=false;
            for (dcount=0; dcount<fdim-1; dcount++) {
              i3=(uvadd) i4/dcountdim1[dcount];
              arSet1[dcount].indx=i3;
              //printf("i3 %d l2 %d\n",i3,l2);
              i4=i4-i3*dcountdim1[dcount];
              if(varsize<=fdim-1) {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[arSet1[dcount].setid].begadd+i3].setsh[varsupsetid[i1]]*varantidim[i1];
                      //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                    break;
                  }
                }
              } else {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[arSet1[dcount].setid].begadd+i3].setsh[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                  }
                }
              }
            }
      /*for (dcount=0; dcount<fdim-1; dcount++) {
        i3=(uvadd) i4/dcountdim1[dcount];
        arSet[dcount].indx=i3;
        i4=i4-i3*dcountdim1[dcount];
        for(i1=0; i1<varsize; i1++) {
          if(vararset[i1]-1==dcount) {
            if(varsubset[i1]==1) {
              l2=l2+ha_setele[ha_set[arSet[dcount].setid].begadd+i3].setsh[varsupsetid[i1]]*varantidim[i1];
            } else {
              l2=l2+i3*varantidim[i1];
            }
            break;
          }
        }
      }*/
      //printf("varval %f i %d\n",ha_cofvar[begadd+l2].varval,begadd+l2);
      temp2=ha_cofvar[begadd+l2].varval;
      temp1=ha_newfpcal(ha_cofvar,ha_set,ha_setele,ha_sumele,ha_calvar1,ha_calvarsize,arSet1,fdim-1,zerodivide);
      if(temp1-ha_cofvar[begadd+l2].varval>0.000000001||temp1-ha_cofvar[begadd+l2].varval<-0.000000001)ha_cofvar[begadd+l2].varval=ha_cofvar[begadd+l2].csolpupd+2*(temp1-ha_cofvar[begadd+l2].varval);
      ha_cofvar[begadd+l2].csolpupd=temp2;
      //printf("varval %f i %d\n",ha_cofvar[begadd+l2].varval,begadd+l2);
    }
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
    free(sum_cof);
    free(ha_sumele);
    free(arSet);
    free(ha_calvar);

        if(ha_cof[index].gltype>0){
        #pragma omp parallel private(l) shared(ha_cofvar,ha_cof,ncof,begadd,index)
        {
          if(ha_cof[index].gltype==1){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval<ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 1!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
          if(ha_cof[index].gltype==2){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval<=ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 2!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
          if(ha_cof[index].gltype==3){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval>ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 3!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
          if(ha_cof[index].gltype==4){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval>=ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 4!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
        }
        }
    
  }
  fclose(filehandle);
  return j;
}

uvadd hnew_graggupd(char *fname,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele) {
  FILE * filehandle;
  char commsyntax[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],linecopy[TABREADLINE];
  char vname[NAMESIZE],sumsyntax[NAMESIZE],argu[NAMESIZE];
  char *readitem=NULL,*p=NULL;
  uvadd i,i1,i3,i4,l,l2=0,j=0,nsumele,dcountdim1[4*MAXVARDIM],nloops;//m,
  uvdim fdim,dcount,sup,varsupsetid[MAXVARDIM];
  int ha_calvarsize=0,totalsum,sumcount=1,npow,nmul,ndiv,nplu,nmin,npar,sumindx;
  bool IsChange=false,IsExplicit=false;
  ha_cgetype zerodivide=0,temp1,temp2;
  uvadd varantidim[MAXVARDIM],varsubset[MAXVARDIM],vararset[MAXVARDIM];
  //printf("line %s cofvar %f\n",line,ha_cofvar[23].varval);
  ha_cgesetindx *arSet1=NULL;
  hcge_calvars *ha_calvar1= NULL;
  strcpy(commsyntax,"update");
  filehandle = fopen(fname,"r");
  while (ha_cgertabl1(commsyntax,filehandle,line,ha_cofvar,ha_cof,ncof,&zerodivide,TABREADLINE)) {
    IsChange=false;
    IsExplicit=false;
    if(strstr(line, "(change)")!=NULL) {
      IsChange=true;
      ha_cgefrstr1(line, "(change)", "");
    }
    if(strstr(line, "(explicit)")!=NULL) {
      IsExplicit=true;
      ha_cgefrstr1(line, "(explicit)", "");
    }
    ha_cgefrstr1(line, commsyntax, "");
    while (ha_cgefrstr(line," ", ""));
    while (ha_cgefrchr(line, '[', '('));
    while (ha_cgefrchr(line, ']', ')'));
    while (ha_cgefrchr(line, '{', '('));
    while (ha_cgefrchr(line, '}', ')'));
    strcpy(linecopy,line);
    //printf("line %s cofvar %f\n",line,ha_cofvar[23].varval);
    //printf("line %s\n",line);
    readitem = strtok(line,"=");
    fdim=ha_cgenchf(readitem, '(');
    if (fdim==1) {
      fdim=fdim+1;
    }
    ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (fdim+1,sizeof(ha_cgesetindx));

    nloops=1;
    if (fdim==0) {
      readitem = strtok(line,"=");
    } else {
      for (i=0; i<fdim-1; i++) {
        if(i==0) {
          readitem = strtok(line,",");
        } else {
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
            //} else {
            //  arSet[i].SuperSetSize=ha_set[i4].size;
            //}
            break;
          }
        nloops=nloops*ha_set[arSet[i].setid].size;
      }
      readitem = strtok(NULL,"=");
      dcountdim1[fdim-2]=1;
      for (i=fdim-3; i>-1; i--) {
        dcountdim1[i]=ha_set[arSet[i+1].setid].size*dcountdim1[i+1];
      }
    }
    strcpy(vname,readitem);
    //printf("vname %s\n",vname);
    strcpy(line,linecopy);
    readitem = strtok(line,"=");
    readitem = strtok(NULL,";");
    if(IsChange==false&&IsExplicit==false) {
      //if(lastupd){
      //while (ha_cgefrstr(readitem,"*", "/100)!(1+"));
      //while (ha_cgefrstr(readitem,"!", "*"));
      //strcpy(line1,vname);
      //strcat(line1,"*(1+");
      //strcat(line1,readitem);
      //strcat(line1,"/100)");
      //readitem=line1;
      //} else{
      while (ha_cgefrstr(readitem,"*", "+"));
      strcpy(line1,vname);
      strcat(line1,"*(1+(");
      strcat(line1,readitem);
      strcat(line1,")/100)");
      readitem=line1;
      //}
    }
    if(IsChange==true) {
      strcpy(line1,vname);
      if(readitem[0]!='+'||readitem[0]!='-') strcat(line1,"+");
      strcat(line1,readitem);
      readitem=line1;
    }
    npow=ha_cgenchf(readitem, '^');
    nmul=ha_cgenchf(readitem, '*');
    ndiv=ha_cgenchf(readitem, '/');
    nmul=nmul+ndiv;
    nplu=ha_cgenchf(readitem, '+');
    nmin=ha_cgenchf(readitem, '-');
    nplu=nplu+nmin;

    //strcpy(line,linecopy);
    //readitem = strtok(line,"=");
    //readitem = strtok(NULL,";");
    while (ha_cgerecovar(readitem)==1);
    //printf("line1 %s\n",readitem);
    hnew_intrpl(readitem);
    //printf("line2 %s\n",readitem);
    npar=ha_cgenchf(readitem, '(');
    strcpy(sumsyntax,"sum(");
    totalsum=hcge_nsum(readitem,sumsyntax);
    hcge_sumcof *sum_cof= (hcge_sumcof *) calloc (totalsum,sizeof(hcge_sumcof));
    sumcount=0;
    strcpy(line1,readitem);
    //printf("line1 %s\n",line1);
    strcpy(line2,line1);
    readitem=line2;
    while (hcge_dsum(readitem,sumsyntax,sum_cof,arSet,ha_set,nset,fdim,sumcount)==1) {
      sumcount++;
    }
    totalsum=sumcount;
    i3=0;
    for (i=0; i<totalsum; i++) {
      i1=1;
      for(j=0; j<sum_cof[i].size; j++) {
        i1=i1*ha_set[sum_cof[i].setid[j]].size;
      }
      sum_cof[i].begadd=i3;
      i3=i3+i1;
    }
    nsumele=i3;
    for (i=0; i<totalsum; i++) {
      i1=1;
      sum_cof[i].antidims[sum_cof[i].size-1]=1;
      for(j=sum_cof[i].size-2; j>-1; j--) {
        sum_cof[i].antidims[j]=sum_cof[i].antidims[j+1]*ha_set[sum_cof[i].setid[j+1]].size;
      }
    }
    hcge_calvars *ha_calvar= (hcge_calvars *) calloc (npow+nmul+nplu+2*npar+2,sizeof(hcge_calvars));
    ha_cgesumele *ha_sumele= (ha_cgesumele *) calloc (nsumele,sizeof(ha_cgesumele));
    //printf("read2 %s\n",readitem);
    //hcge_rsumele(sum_cof,totalsum,ha_set,nset,ha_setele,ha_sumele);
    //printf("read1 %s\n",readitem);
    sumcount=0;
    strcpy(line2,line1);
    readitem=line2;
    sumindx=0;
    while (hnew_calsum(readitem,sumsyntax,ha_set,nset,ha_setele,ha_cofvar,ncofvar,ncofele,ha_cof,ncof,ha_var,nvar,sum_cof,totalsum,ha_sumele,nsumele,ha_calvar,arSet,fdim,&sumindx,sumcount,zerodivide)==1) {
      sumcount++;
    }
    //printf("read %s\n",readitem);
    //for (i=0;i<nsumele;i++) printf("varele %s varval%f\n",ha_sumele[i].varname,ha_sumele[i].varval);
    strcpy(line1,readitem);
    uvadd index=ncof-1, begadd=0;//,simpl=0;
    bool check10=true;
    //if (strpbrk(readitem,"*+-^/)")==NULL) {
    //simpl=1;
    //}
        uvadd varsize=0;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            begadd=ha_cof[index].begadd;
            varsize=ha_cof[index].size;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            check10=false;
            break;
          }
        } while (index--);
        if (check10) {
          index=nvar-1;
          do {
            if (strcmp(ha_var[index].cofname,p)==0) {
              begadd=ncofele+ha_var[index].begadd;
              varsize=ha_var[index].size;
              if(ha_var[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
              break;
            }
          } while (index--);
        }
        //if(strcmp("rgdpcons",p)==0){
        //printf("line %s\n",linecopy);
        //for (i=0;i<totalsum;i++) printf("name %s\n",sum_cof[i].sumname);
        //for (i=0;i<totalsum;i++)for(j=0;j<sum_cof[i].summatsize;j++) printf("varname %s varval %f\n",sum_cof[i].sumname,ha_sumele[sum_cof[i].begadd+j].varval);
        //for (i=0;i<totalsum;i++){printf("name %s ",sum_cof[i].sumname);for(j=0;j<sum_cof[i].size;j++)printf("set %s ",ha_set[sum_cof[i].setid[j]].setname);printf("\n");}
        //}
        //uvadd *varantidim= (uvadd *) calloc (varsize,sizeof(uvadd));
        //uvadd *varsubset= (uvadd *) calloc (varsize,sizeof(uvadd));
        //uvadd *vararset= (uvadd *) calloc (varsize,sizeof(uvadd));
        //printf("argu OK!!!! %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){varantidim[l]=0;varsubset[l]=0;varsupsetid[l]=0;}
        if (check10) {
            for (dcount=0; dcount<ha_var[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  varantidim[dcount]=ha_var[index].antidims[dcount];
                  vararset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_var[index].setid[dcount]].size){
                    //varsubset[dcount]=0;
                  //} else {
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                    //printf("arset %s varsupsetid %d\n",ha_set[arSet[l].setid].setname,varsupsetid[dcount]);
                  }
                  break;
                }
            }
        } else {
            //printf("cof %s argu %s\n",ha_cof[index].cofname,argu);
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  varantidim[dcount]=ha_cof[index].antidims[dcount];
                  vararset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    //varsubset[dcount]=0;
                  //} else {
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                    //printf("arset %s sup0 %d sup1 %s sup2 %d\n",ha_set[arSet[l].setid].setname,ha_set[arSet[l].setid].subsetid[0],ha_set[ha_set[arSet[l].setid].subsetid[1]].setname,ha_set[arSet[l].setid].subsetid[2]);
                    //printf("arset %s supset %d varsupsetid %d cof set size %d\n",ha_set[arSet[l].setid].setname,ha_set[arSet[l].setid].subsetid[sup],varsupsetid[dcount],ha_cof[index].setid[dcount]);
                    //printf("sup %d\n",sup);
                  }
                  break;
                }
            }
        }
    /*uvadd varsize=0;
    p=strtok(vname,"(");
    do {
      if (strcmp(ha_cof[index].cofname,p)==0) {
        begadd=ha_cof[index].begadd;
        varsize=ha_cof[index].size;
        check10=false;
        break;
      }
    } while (index--);
    if (check10) {
      index=nvar-1;
      do {
        if (strcmp(ha_var[index].cofname,p)==0) {
          begadd=ncofele+ha_var[index].begadd;
          varsize=ha_var[index].size;
          break;
        }
      } while (index--);
    }
    uvadd *varantidim= (uvadd *) calloc (varsize,sizeof(uvadd));
    uvadd *varsubset= (uvadd *) calloc (varsize,sizeof(uvadd));
    uvadd *vararset= (uvadd *) calloc (varsize,sizeof(uvadd));
    if (check10&&ha_var[index].size>0) {
      for (dcount=0; dcount<ha_var[index].size-1; dcount++) {
        p=strtok(NULL,",");
        for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
            varantidim[dcount]=ha_var[index].antidims[dcount];
            vararset[dcount]=l+1;
            if (strcmp(ha_set[arSet[l].setid].setname,ha_set[ha_var[index].setid[dcount]].setname)==0) {
              varsubset[dcount]=0;
            } else {
              varsubset[dcount]=1;
              for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupid[dcount]=sup;break;}
            }
            break;
          }
      }
      p=strtok(NULL,")");
      for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
          varantidim[dcount]=ha_var[index].antidims[dcount];
          vararset[dcount]=l+1;
          if (strcmp(ha_set[arSet[l].setid].setname,ha_set[ha_var[index].setid[dcount]].setname)==0) {
            varsubset[dcount]=0;
          } else {
            varsubset[dcount]=1;
            for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupid[dcount]=sup;break;}
          }
          break;
        }
    }
    if (!check10&&ha_cof[index].size>0) {
      for (dcount=0; dcount<ha_cof[index].size-1; dcount++) {
        p=strtok(NULL,",");
        for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
            varantidim[dcount]=ha_cof[index].antidims[dcount];
            vararset[dcount]=l+1;
            if (strcmp(ha_set[arSet[l].setid].setname,ha_set[ha_cof[index].setid[dcount]].setname)==0) {
              varsubset[dcount]=0;
            } else {
              varsubset[dcount]=1;
              for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupid[dcount]=sup;break;}
            }
            break;
          }
      }
      p=strtok(NULL,")");
      for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
          varantidim[dcount]=ha_cof[index].antidims[dcount];
          vararset[dcount]=l+1;
          if (strcmp(ha_set[arSet[l].setid].setname,ha_set[ha_cof[index].setid[dcount]].setname)==0) {
            varsubset[dcount]=0;
          } else {
            varsubset[dcount]=1;
            for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupid[dcount]=sup;break;}
          }
          break;
        }
    }*/
    //printf("read %s\n line1 %s fdim %d\n\n",readitem,line1,fdim);
    ha_newfparse(line1,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,&ha_calvarsize,arSet,fdim-1);
    //for(l=0; l<ha_calvarsize; l++) for(l2=0; l2<10; l2++) {
    //printf("name %s oper %d tmptype1 %d begadd1 %d ad1 %d tmptype2 %d begadd2 %d ad2 %d\n",ha_calvar[l].TmpVarName,ha_calvar[l].Oper,ha_calvar[l].Var1Type,ha_calvar[l].Var1BegAdd,ha_calvar[l].Var1ADims[l2],ha_calvar[l].Var2Type,ha_calvar[l].Var2BegAdd,ha_calvar[l].Var2ADims[l2]);
    //}
        #pragma omp parallel private(l,l2,i4,dcount,i3,i1,temp1,temp2,arSet1,ha_calvar1) shared(ha_cofvar,arSet)
        {
        if(omp_get_thread_num()!=0){
          arSet1=realloc(arSet1,(fdim+1)*sizeof(ha_cgesetindx));
          memcpy (arSet1,arSet,(fdim+1)*sizeof(ha_cgesetindx));
          ha_calvar1=realloc(ha_calvar1,ha_calvarsize*sizeof(hcge_calvars));
          memcpy (ha_calvar1,ha_calvar,ha_calvarsize*sizeof(hcge_calvars));
        }else{
          ha_calvar1=ha_calvar;
          arSet1=arSet;
        }
        #pragma omp for
    for (l=0; l<nloops; l++) {
      l2=0;
      i4=l;
      //check=false;
            for (dcount=0; dcount<fdim-1; dcount++) {
              i3=(uvadd) i4/dcountdim1[dcount];
              arSet1[dcount].indx=i3;
              //printf("i3 %d l2 %d\n",i3,l2);
              i4=i4-i3*dcountdim1[dcount];
              if(varsize<=fdim-1) {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[arSet1[dcount].setid].begadd+i3].setsh[varsupsetid[i1]]*varantidim[i1];
                      //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                    break;
                  }
                }
              } else {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[arSet1[dcount].setid].begadd+i3].setsh[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                  }
                }
              }
            }
      /*for (dcount=0; dcount<fdim-1; dcount++) {
        i3=(uvadd) i4/dcountdim1[dcount];
        arSet[dcount].indx=i3;
        i4=i4-i3*dcountdim1[dcount];
        for(i1=0; i1<varsize; i1++) {
          if(vararset[i1]-1==dcount) {
            if(varsubset[i1]==1) {
              l2=l2+ha_setele[ha_set[arSet[dcount].setid].begadd+i3].setsh[varsupid[i1]]*varantidim[i1];
            } else {
              l2=l2+i3*varantidim[i1];
            }
            break;
          }
        }
      }*/
      //printf("varval %f i %d\n",ha_cofvar[begadd+l2].varval,begadd+l2);
      //temp2=ha_cofvar[begadd+l2].varval;
      temp1=ha_newfpcal(ha_cofvar,ha_set,ha_setele,ha_sumele,ha_calvar1,ha_calvarsize,arSet1,fdim-1,zerodivide);
      if(temp1-ha_cofvar[begadd+l2].varval>0.000000001||temp1-ha_cofvar[begadd+l2].varval<-0.000000001)ha_cofvar[begadd+l2].varval=0.5*(ha_cofvar[begadd+l2].csolpupd+temp1);
      //ha_cofvar[begadd+l2].csolpupd=temp2;
      //printf("varval %f i %d\n",ha_cofvar[begadd+l2].varval,begadd+l2);
    }
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
    free(sum_cof);
    free(ha_sumele);
    //free(varantidim);
    //free(varsubset);
    free(arSet);
    free(ha_calvar);

        if(ha_cof[index].gltype>0){
        #pragma omp parallel private(l) shared(ha_cofvar,ha_cof,ncof,begadd,index)
        {
          if(ha_cof[index].gltype==1){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval<ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 1!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
          if(ha_cof[index].gltype==2){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval<=ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 2!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
          if(ha_cof[index].gltype==3){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval>ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 3!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
          if(ha_cof[index].gltype==4){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval>=ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 4!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
        }
        }
    
  }
  fclose(filehandle);
  return j;
}

uvadd hnew_gupd(char *fname,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele) {
  FILE * filehandle;
  char commsyntax[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],linecopy[TABREADLINE];
  char vname[NAMESIZE],sumsyntax[NAMESIZE],argu[NAMESIZE];
  char *readitem=NULL,*p=NULL;
  uvadd i,i1,i3,i4,l,l2=0,j=0,nsumele,dcountdim1[4*MAXVARDIM],nloops;//m,
  uvdim fdim,dcount,sup,varsupsetid[MAXVARDIM];
  int ha_calvarsize=0,totalsum,sumcount=1,npow,nmul,ndiv,nplu,nmin,npar,sumindx;
  bool IsChange=false,IsExplicit=false;
  ha_cgetype zerodivide=0,temp1;
  uvadd varantidim[MAXVARDIM],varsubset[MAXVARDIM],vararset[MAXVARDIM];
  //printf("line %s cofvar %f\n",line,ha_cofvar[23].varval);
  ha_cgesetindx *arSet1=NULL;
  hcge_calvars *ha_calvar1= NULL;
  strcpy(commsyntax,"update");
  filehandle = fopen(fname,"r");
  while (ha_cgertabl1(commsyntax,filehandle,line,ha_cofvar,ha_cof,ncof,&zerodivide,TABREADLINE)) {
    IsChange=false;
    IsExplicit=false;
    if(strstr(line, "(change)")!=NULL) {
      IsChange=true;
      ha_cgefrstr1(line, "(change)", "");
    }
    if(strstr(line, "(explicit)")!=NULL) {
      IsExplicit=true;
      ha_cgefrstr1(line, "(explicit)", "");
    }
    ha_cgefrstr1(line, commsyntax, "");
    while (ha_cgefrstr(line," ", ""));
    while (ha_cgefrchr(line, '[', '('));
    while (ha_cgefrchr(line, ']', ')'));
    while (ha_cgefrchr(line, '{', '('));
    while (ha_cgefrchr(line, '}', ')'));
    strcpy(linecopy,line);
    //printf("line %s cofvar %f\n",line,ha_cofvar[23].varval);
    //printf("line %s\n",line);
    readitem = strtok(line,"=");
    fdim=ha_cgenchf(readitem, '(');
    if (fdim==1) {
      fdim=fdim+1;
    }
    ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (fdim+1,sizeof(ha_cgesetindx));

    nloops=1;
    if (fdim==0) {
      readitem = strtok(line,"=");
    } else {
      for (i=0; i<fdim-1; i++) {
        if(i==0) {
          readitem = strtok(line,",");
        } else {
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
            //} else {
            //  arSet[i].SuperSetSize=ha_set[i4].size;
            //}
            break;
          }
        nloops=nloops*ha_set[arSet[i].setid].size;
      }
      readitem = strtok(NULL,"=");
      dcountdim1[fdim-2]=1;
      for (i=fdim-3; i>-1; i--) {
        dcountdim1[i]=ha_set[arSet[i+1].setid].size*dcountdim1[i+1];
      }
    }
    strcpy(vname,readitem);
    //printf("vname %s\n",vname);
    strcpy(line,linecopy);
    readitem = strtok(line,"=");
    readitem = strtok(NULL,";");
    if(IsChange==false&&IsExplicit==false) {
      //if(lastupd){
      while (ha_cgefrstr(readitem,"*", "/100)!(1+"));
      while (ha_cgefrstr(readitem,"!", "*"));
      strcpy(line1,vname);
      strcat(line1,"*(1+");
      strcat(line1,readitem);
      strcat(line1,"/100)");
      readitem=line1;
      //} else{
      //while (ha_cgefrstr(readitem,"*", "+"));
      //strcpy(line1,vname);
      //strcat(line1,"*(1+(");
      //strcat(line1,readitem);
      //strcat(line1,")/100)");
      //readitem=line1;
      //}
    }
    if(IsChange==true) {
      strcpy(line1,vname);
      if(readitem[0]!='+'||readitem[0]!='-') strcat(line1,"+");
      strcat(line1,readitem);
      readitem=line1;
    }
    npow=ha_cgenchf(readitem, '^');
    nmul=ha_cgenchf(readitem, '*');
    ndiv=ha_cgenchf(readitem, '/');
    nmul=nmul+ndiv;
    nplu=ha_cgenchf(readitem, '+');
    nmin=ha_cgenchf(readitem, '-');
    nplu=nplu+nmin;

    //strcpy(line,linecopy);
    //readitem = strtok(line,"=");
    //readitem = strtok(NULL,";");
    while (ha_cgerecovar(readitem)==1);
    //printf("line1 %s\n",readitem);
    hnew_intrpl(readitem);
    //printf("line2 %s\n",readitem);
    npar=ha_cgenchf(readitem, '(');
    strcpy(sumsyntax,"sum(");
    totalsum=hcge_nsum(readitem,sumsyntax);
    hcge_sumcof *sum_cof= (hcge_sumcof *) calloc (totalsum,sizeof(hcge_sumcof));
    sumcount=0;
    strcpy(line1,readitem);
    //printf("line1 %s\n",line1);
    strcpy(line2,line1);
    readitem=line2;
    while (hcge_dsum(readitem,sumsyntax,sum_cof,arSet,ha_set,nset,fdim,sumcount)==1) {
      sumcount++;
    }
    totalsum=sumcount;
    i3=0;
    for (i=0; i<totalsum; i++) {
      i1=1;
      for(j=0; j<sum_cof[i].size; j++) {
        i1=i1*ha_set[sum_cof[i].setid[j]].size;
      }
      sum_cof[i].begadd=i3;
      i3=i3+i1;
    }
    nsumele=i3;
    for (i=0; i<totalsum; i++) {
      i1=1;
      sum_cof[i].antidims[sum_cof[i].size-1]=1;
      for(j=sum_cof[i].size-2; j>-1; j--) {
        sum_cof[i].antidims[j]=sum_cof[i].antidims[j+1]*ha_set[sum_cof[i].setid[j+1]].size;
      }
    }
    hcge_calvars *ha_calvar= (hcge_calvars *) calloc (npow+nmul+nplu+2*npar+2,sizeof(hcge_calvars));
    ha_cgesumele *ha_sumele= (ha_cgesumele *) calloc (nsumele,sizeof(ha_cgesumele));
    //printf("read2 %s\n",readitem);
    //hcge_rsumele(sum_cof,totalsum,ha_set,nset,ha_setele,ha_sumele);
    //printf("read1 %s\n",readitem);
    sumcount=0;
    strcpy(line2,line1);
    readitem=line2;
    sumindx=0;
    while (hnew_calsum(readitem,sumsyntax,ha_set,nset,ha_setele,ha_cofvar,ncofvar,ncofele,ha_cof,ncof,ha_var,nvar,sum_cof,totalsum,ha_sumele,nsumele,ha_calvar,arSet,fdim,&sumindx,sumcount,zerodivide)==1) {
      sumcount++;
    }
    //printf("read %s\n",readitem);
    //for (i=0;i<nsumele;i++) printf("varele %s varval%f\n",ha_sumele[i].varname,ha_sumele[i].varval);
    strcpy(line1,readitem);
    uvadd index=ncof-1, begadd=0;//,simpl=0;
    bool check10=true;
    //if (strpbrk(readitem,"*+-^/)")==NULL) {
    //simpl=1;
    //}
        uvadd varsize=0;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            begadd=ha_cof[index].begadd;
            varsize=ha_cof[index].size;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            check10=false;
            break;
          }
        } while (index--);
        if (check10) {
          index=nvar-1;
          do {
            if (strcmp(ha_var[index].cofname,p)==0) {
              begadd=ncofele+ha_var[index].begadd;
              varsize=ha_var[index].size;
              if(ha_var[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
              break;
            }
          } while (index--);
        }
        //if(strcmp("rgdpcons",p)==0){
        //printf("line %s\n",linecopy);
        //for (i=0;i<totalsum;i++) printf("name %s\n",sum_cof[i].sumname);
        //for (i=0;i<totalsum;i++)for(j=0;j<sum_cof[i].summatsize;j++) printf("varname %s varval %f\n",sum_cof[i].sumname,ha_sumele[sum_cof[i].begadd+j].varval);
        //for (i=0;i<totalsum;i++){printf("name %s ",sum_cof[i].sumname);for(j=0;j<sum_cof[i].size;j++)printf("set %s ",ha_set[sum_cof[i].setid[j]].setname);printf("\n");}
        //}
        //uvadd *varantidim= (uvadd *) calloc (varsize,sizeof(uvadd));
        //uvadd *varsubset= (uvadd *) calloc (varsize,sizeof(uvadd));
        //uvadd *vararset= (uvadd *) calloc (varsize,sizeof(uvadd));
        //printf("argu OK!!!! %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){varantidim[l]=0;varsubset[l]=0;varsupsetid[l]=0;}
        if (check10) {
            for (dcount=0; dcount<ha_var[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  varantidim[dcount]=ha_var[index].antidims[dcount];
                  vararset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_var[index].setid[dcount]].size){
                    //varsubset[dcount]=0;
                  //} else {
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                    //printf("arset %s varsupsetid %d\n",ha_set[arSet[l].setid].setname,varsupsetid[dcount]);
                  }
                  break;
                }
            }
        } else {
            //printf("cof %s argu %s\n",ha_cof[index].cofname,argu);
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  varantidim[dcount]=ha_cof[index].antidims[dcount];
                  vararset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    //varsubset[dcount]=0;
                  //} else {
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                    //printf("arset %s sup0 %d sup1 %s sup2 %d\n",ha_set[arSet[l].setid].setname,ha_set[arSet[l].setid].subsetid[0],ha_set[ha_set[arSet[l].setid].subsetid[1]].setname,ha_set[arSet[l].setid].subsetid[2]);
                    //printf("arset %s supset %d varsupsetid %d cof set size %d\n",ha_set[arSet[l].setid].setname,ha_set[arSet[l].setid].subsetid[sup],varsupsetid[dcount],ha_cof[index].setid[dcount]);
                    //printf("sup %d\n",sup);
                  }
                  break;
                }
            }
        }
    /*uvadd varsize=0;
    p=strtok(vname,"(");
    do {
      if (strcmp(ha_cof[index].cofname,p)==0) {
        begadd=ha_cof[index].begadd;
        varsize=ha_cof[index].size;
        check10=false;
        break;
      }
    } while (index--);
    if (check10) {
      index=nvar-1;
      do {
        if (strcmp(ha_var[index].cofname,p)==0) {
          begadd=ncofele+ha_var[index].begadd;
          varsize=ha_var[index].size;
          break;
        }
      } while (index--);
    }
    uvadd *varantidim= (uvadd *) calloc (varsize,sizeof(uvadd));
    uvadd *varsubset= (uvadd *) calloc (varsize,sizeof(uvadd));
    uvadd *vararset= (uvadd *) calloc (varsize,sizeof(uvadd));
    if (check10&&ha_var[index].size>0) {
      for (dcount=0; dcount<ha_var[index].size-1; dcount++) {
        p=strtok(NULL,",");
        for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
            varantidim[dcount]=ha_var[index].antidims[dcount];
            vararset[dcount]=l+1;
            if (strcmp(ha_set[arSet[l].setid].setname,ha_set[ha_var[index].setid[dcount]].setname)==0) {
              varsubset[dcount]=0;
            } else {
              varsubset[dcount]=1;
              for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupid[dcount]=sup;break;}
            }
            break;
          }
      }
      p=strtok(NULL,")");
      for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
          varantidim[dcount]=ha_var[index].antidims[dcount];
          vararset[dcount]=l+1;
          if (strcmp(ha_set[arSet[l].setid].setname,ha_set[ha_var[index].setid[dcount]].setname)==0) {
            varsubset[dcount]=0;
          } else {
            varsubset[dcount]=1;
            for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupid[dcount]=sup;break;}
          }
          break;
        }
    }
    if (!check10&&ha_cof[index].size>0) {
      for (dcount=0; dcount<ha_cof[index].size-1; dcount++) {
        p=strtok(NULL,",");
        for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
            varantidim[dcount]=ha_cof[index].antidims[dcount];
            vararset[dcount]=l+1;
            if (strcmp(ha_set[arSet[l].setid].setname,ha_set[ha_cof[index].setid[dcount]].setname)==0) {
              varsubset[dcount]=0;
            } else {
              varsubset[dcount]=1;
              for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupid[dcount]=sup;break;}
            }
            break;
          }
      }
      p=strtok(NULL,")");
      for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
          varantidim[dcount]=ha_cof[index].antidims[dcount];
          vararset[dcount]=l+1;
          if (strcmp(ha_set[arSet[l].setid].setname,ha_set[ha_cof[index].setid[dcount]].setname)==0) {
            varsubset[dcount]=0;
          } else {
            varsubset[dcount]=1;
            for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){varsupid[dcount]=sup;break;}
          }
          break;
        }
    }*/
    //printf("read %s\n line1 %s fdim %d\n\n",readitem,line1,fdim);
    ha_newfparse(line1,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,&ha_calvarsize,arSet,fdim-1);
    //for(l=0; l<ha_calvarsize; l++) for(l2=0; l2<10; l2++) {
    //printf("name %s oper %d tmptype1 %d begadd1 %d ad1 %d tmptype2 %d begadd2 %d ad2 %d\n",ha_calvar[l].TmpVarName,ha_calvar[l].Oper,ha_calvar[l].Var1Type,ha_calvar[l].Var1BegAdd,ha_calvar[l].Var1ADims[l2],ha_calvar[l].Var2Type,ha_calvar[l].Var2BegAdd,ha_calvar[l].Var2ADims[l2]);
    //}
        #pragma omp parallel private(l,l2,i4,dcount,i3,i1,temp1,arSet1,ha_calvar1) shared(ha_cofvar,arSet)
        {
        if(omp_get_thread_num()!=0){
          arSet1=realloc(arSet1,(fdim+1)*sizeof(ha_cgesetindx));
          memcpy (arSet1,arSet,(fdim+1)*sizeof(ha_cgesetindx));
          ha_calvar1=realloc(ha_calvar1,ha_calvarsize*sizeof(hcge_calvars));
          memcpy (ha_calvar1,ha_calvar,ha_calvarsize*sizeof(hcge_calvars));
        }else{
          ha_calvar1=ha_calvar;
          arSet1=arSet;
        }
        #pragma omp for
    for (l=0; l<nloops; l++) {
      l2=0;
      i4=l;
      //check=false;
            for (dcount=0; dcount<fdim-1; dcount++) {
              i3=(uvadd) i4/dcountdim1[dcount];
              arSet1[dcount].indx=i3;
              //printf("i3 %d l2 %d\n",i3,l2);
              i4=i4-i3*dcountdim1[dcount];
              if(varsize<=fdim-1) {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[arSet1[dcount].setid].begadd+i3].setsh[varsupsetid[i1]]*varantidim[i1];
                      //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                    break;
                  }
                }
              } else {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[arSet1[dcount].setid].begadd+i3].setsh[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                  }
                }
              }
            }
      /*for (dcount=0; dcount<fdim-1; dcount++) {
        i3=(uvadd) i4/dcountdim1[dcount];
        arSet[dcount].indx=i3;
        i4=i4-i3*dcountdim1[dcount];
        for(i1=0; i1<varsize; i1++) {
          if(vararset[i1]-1==dcount) {
            if(varsubset[i1]==1) {
              l2=l2+ha_setele[ha_set[arSet[dcount].setid].begadd+i3].setsh[varsupid[i1]]*varantidim[i1];
            } else {
              l2=l2+i3*varantidim[i1];
            }
            break;
          }
        }
      }*/
      //printf("varval %f i %d\n",ha_cofvar[begadd+l2].varval,begadd+l2);
      temp1=ha_newfpcal(ha_cofvar,ha_set,ha_setele,ha_sumele,ha_calvar1,ha_calvarsize,arSet1,fdim-1,zerodivide);
      if(temp1-ha_cofvar[begadd+l2].varval>0.000000001||temp1-ha_cofvar[begadd+l2].varval<-0.000000001)ha_cofvar[begadd+l2].varval=temp1;
      //printf("varval %f i %d\n",ha_cofvar[begadd+l2].varval,begadd+l2);
    }
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
    free(sum_cof);
    free(ha_sumele);
    //free(varantidim);
    //free(varsubset);
    free(arSet);
    free(ha_calvar);

        if(ha_cof[index].gltype>0){
        #pragma omp parallel private(l) shared(ha_cofvar,ha_cof,ncof,begadd,index)
        {
          if(ha_cof[index].gltype==1){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval<ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 1!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
          if(ha_cof[index].gltype==2){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval<=ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 2!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
          if(ha_cof[index].gltype==3){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval>ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 3!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
          if(ha_cof[index].gltype==4){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(ha_cofvar[begadd+l].varval>=ha_cof[index].glval){
              printf("Error!!! Condition not met var %s type 4!\n",ha_cof[index].cofname);
              l=varsize;
            }
          }
          }
        }
        }
    
  }
  fclose(filehandle);
  return j;
}

int hnew_calsum(char *formulain, char *commsyntax,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele,ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele, hcge_cof *ha_cof,uvadd ncof, hcge_cof *ha_var,uvadd nvar,hcge_sumcof *sum_cof,int totalsum,ha_cgesumele *ha_sumele,uvadd nsumele,hcge_calvars *ha_calvar,ha_cgesetindx *arSet1,uvdim fdim,int *sumindx,int j, ha_cgetype zerodivide) {
  char *readitem,*p;//,*p1,interchar2[NAMESIZE],line5[TABREADLINE];
  char interchar[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE];//,line3[TABREADLINE],line4[TABREADLINE];//,interchar1[NAMESIZE]
  int ha_calvarsize,length,k=0,k1=0,i=0;
  uvdim fdimsumcof,setsh,dcount;
  uvadd l,l1,l2,nloops,dcountdim1[4*MAXVARDIM];
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
      k1=ha_cgefind(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      } else {
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
          } else {
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
        //printf("p %s hacovar %f\n",p,ha_cofvar[23].varval);
        ha_newfparse(p,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,&ha_calvarsize,arSet,fdimsumcof);
        //for(l=0; l<ha_calvarsize; l++) for(l2=0; l2<10; l2++) {
        //printf("name %s oper %d tmptype %d begadd %d ad %d\n",ha_calvar[l].TmpVarName,ha_calvar[l].Oper,ha_calvar[l].Var1Type,ha_calvar[l].Var1BegAdd,ha_calvar[l].Var1ADims[l2]);
        //}
        //for(l=0;l<ha_calvarsize;l++) for(l2=0;l2<10;l2++) printf("tmptype %d begadd %d ad %d\n",ha_calvar[l].TmpVarType,ha_calvar[l].Var1BegAdd,ha_calvar[l].Var1ADims[l2]);
        #pragma omp parallel private(l,l1,l2,dcount,setsh,vval,arSet2,ha_calvar1) shared(ha_cofvar,arSet,ha_sumele)
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
        for (l=0; l<nloops; l++) {
          l2=l;
          for (dcount=0; dcount<sum_cof[j].size; dcount++) {
            setsh=(uvdim) l2/dcountdim1[dcount];
            arSet2[dcount].indx=setsh;
            l2=l2-setsh*dcountdim1[dcount];
          }
          vval=0;
          for (l1=0; l1<ha_set[sum_cof[j].sumsetid].size; l1++) {
            arSet2[sum_cof[j].size].indx=l1;
            vval+=ha_newfpcal(ha_cofvar,ha_set,ha_setele,ha_sumele,ha_calvar1,ha_calvarsize,arSet2,fdimsumcof,zerodivide);
          }
          ha_sumele[*sumindx+l].varval=vval;
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
        } else {
          if (interchar[strlen(interchar)-1]=='{') {
            interchar[strlen(interchar)-1]='\0';
          } else {
            strcat(interchar,"}");
          }
        }
        while(ha_cgefrstr(formulain,line1,interchar)!=NULL);
        free(arSet);
        return 1;
      }
    } else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      ha_cgecutsum(line);
      k1=ha_cgefind(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      } else {
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
          } else {
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
        //printf("p %s hacovar %f\n",p,ha_cofvar[23].varval);
        ha_newfparse(p,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,&ha_calvarsize,arSet,fdimsumcof);
        //for(l=0;l<sum_cof[j].size+1;l++)printf("name %s set %s\n",sum_cof[j].sumname,ha_set[arSet[l].setid].setname);
        //for(l=0; l<ha_calvarsize; l++) for(l2=0; l2<10; l2++) {
        //printf("name %s oper %d tmptype %d begadd %d ad %d\n",ha_calvar[l].TmpVarName,ha_calvar[l].Oper,ha_calvar[l].Var1Type,ha_calvar[l].Var1BegAdd,ha_calvar[l].Var1ADims[l2]);
        //}
        //for(l=0;l<ha_calvarsize;l++) for(l2=0;l2<10;l2++) printf("tmptype %d begadd %d ad %d\n",ha_calvar[l].TmpVarType,ha_calvar[l].Var1BegAdd,ha_calvar[l].Var1ADims[l2]);
        #pragma omp parallel private(l,l1,l2,dcount,setsh,vval,arSet2,ha_calvar1) shared(ha_cofvar,arSet,ha_sumele)
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
        for (l=0; l<nloops; l++) {
          l2=l;
          for (dcount=0; dcount<sum_cof[j].size; dcount++) {
            setsh=(uvdim) l2/dcountdim1[dcount];
            arSet2[dcount].indx=setsh;
            l2=l2-setsh*dcountdim1[dcount];
          }
          vval=0;
          for (l1=0; l1<ha_set[sum_cof[j].sumsetid].size; l1++) {
            arSet2[sum_cof[j].size].indx=l1;
            //strcpy(p,line2);
            vval+=ha_newfpcal(ha_cofvar,ha_set,ha_setele,ha_sumele,ha_calvar1,ha_calvarsize,arSet2,fdimsumcof,zerodivide);
            //printf("vval %f\n",vval);
            //if (simpl==1) {
            //vval=vval+hcge_varrepl(p,ha_cofvar,ncofele,ha_cof,ncof,ha_var,nvar,ha_setele,sum_cof,totalsum,ha_sumele,nsumele,arSet,fdimsumcof);//hcge_ssimplrpl(p,ha_cofvar,ncofele,ha_cof,ncof,ha_var,nvar,sum_cof,totalsum,ha_sumele,nsumele,arSet,fdimsumcof);
            //} else {
            //vval=vval+ha_cgefparse(p,ha_cofvar,ncofvar,ncofele,ha_cof,ncof,ha_var,nvar,ha_setele,sum_cof,totalsum,ha_sumele,nsumele,vartemppow,vartempmuldiv,vartempplu,arSet,fdimsumcof,zerodivide);
            //}
          }
          ha_sumele[*sumindx+l].varval=vval;//ha_sumele[*sumindx+l2].varval=vval;
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
        } else {
          if (interchar[strlen(interchar)-1]=='{') {
            interchar[strlen(interchar)-1]='\0';
          } else {
            strcat(interchar,"}");
          }
        }
        while(ha_cgefrstr(formulain,line1,interchar));
        free(arSet);
        return 1;
      }
    } else {
      i=i+k+4;
      readitem=formulain+i;
    }
  }
  return 0;
}

uvadd hnew_biupd(PetscInt rank,char *fname,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele,ha_cgeexovar *ha_cgeshock,uvadd nvarele,int laA,uvdim subints,bool IsIni,int IsSplint,int nsteps) {
  FILE * filehandle,*fout;
  char commsyntax[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],linecopy[TABREADLINE];
  char vname[NAMESIZE],sumsyntax[NAMESIZE],argu[NAMESIZE];
  char *readitem=NULL,*p=NULL;
  char filename[1024],j1name[1024];
  uvadd i,i1,i3,i4,l,l2=0,j=0,nsumele,dcountdim1[4*MAXVARDIM],nloops,nloops1,tsize,stosize,bsize,rsize,vvar,xvar,yvar,jvar,wvar,svar,vlmu,dlmu,wlmu,varsize,l2v;//m,
  uvdim fdim,fdim1,dcount,sup;
  bool IsChange=false;
  uvadd index=ncof-1, begadd=0,sbegadd=0,address=0;//,simpl=0;
  bool check10=true;
  ha_cgetype zerodivide=0,temp1,temp2,temp3,maxerr=0,temp4,temp5,sx0,sxn;
  uvadd varantidim[MAXVARDIM],varsubset[MAXVARDIM],vararset[MAXVARDIM],varsupsetid[MAXVARDIM];
  uvadd xantidim[MAXVARDIM],xsubset[MAXVARDIM],xarset[MAXVARDIM],xsupsetid[MAXVARDIM];
  uvadd yantidim[MAXVARDIM],ysubset[MAXVARDIM],yarset[MAXVARDIM],ysupsetid[MAXVARDIM];
  uvadd jantidim[MAXVARDIM],jsubset[MAXVARDIM],jarset[MAXVARDIM],jsupsetid[MAXVARDIM];
  uvadd wantidim[MAXVARDIM],wsubset[MAXVARDIM],warset[MAXVARDIM],wsupsetid[MAXVARDIM];
  uvadd santidim[MAXVARDIM],ssubset[MAXVARDIM],sarset[MAXVARDIM],ssupsetid[MAXVARDIM];
  uvadd mvantidim[MAXVARDIM],mvsubset[MAXVARDIM],mvarset[MAXVARDIM],mvsupsetid[MAXVARDIM];
  uvadd mdantidim[MAXVARDIM],mdsubset[MAXVARDIM],mdarset[MAXVARDIM],mdsupsetid[MAXVARDIM];
  uvadd mwantidim[MAXVARDIM],mwsubset[MAXVARDIM],mwarset[MAXVARDIM],mwsupsetid[MAXVARDIM];
  size_t freadresult;
  //printf("line %s cofvar %f\n",line,ha_cofvar[23].varval);
  ha_cgesetindx *arSet1=NULL;
  hcge_calvars *ha_calvar1= NULL;
  strcpy(commsyntax,"splinter");
      if(rank<10)strcpy(j1name,"000");
      if(rank<100&&rank>9)strcpy(j1name,"00");
      if(rank<1000&&rank>99)strcpy(j1name,"0");
      if(rank>=1000)j1name[0]='\0';
      sprintf(filename, "%d",rank);
      strcat(j1name,filename);
      strcpy(filename,"_biupd");
      strcat(filename,j1name);
      strcat(filename,".bin");
  filehandle = fopen(fname,"r");
  while (ha_cgertabl1(commsyntax,filehandle,line,ha_cofvar,ha_cof,ncof,&zerodivide,TABREADLINE)) {
    IsChange=false;
    if(strstr(line, "(change)")!=NULL) {
      IsChange=true;
      ha_cgefrstr1(line, "(change)", "");
    }
    ha_cgefrstr1(line, commsyntax, "");
    while (ha_cgefrstr(line," ", ""));
    while (ha_cgefrchr(line, '[', '('));
    while (ha_cgefrchr(line, ']', ')'));
    while (ha_cgefrchr(line, '{', '('));
    while (ha_cgefrchr(line, '}', ')'));
    strcpy(linecopy,line);
    //printf("line %s cofvar %f\n",line,ha_cofvar[23].varval);
    //printf("line %s\n",line);
    readitem = strtok(line,"=");//Syntax (linear or change)(other sets):(inter set):(spline set)var_plus()=var():x():x_jump():weight();//
    //printf("read %s\n",readitem);
    fdim=ha_cgenchf(readitem, '(');
    if (fdim==1) {
      fdim=fdim+1;
    }
    fdim1=fdim-1;
    ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (fdim+1,sizeof(ha_cgesetindx));
    //printf("fdim %d\n",fdim);
    if(fdim<3){printf("We should at least have an intertemporal variable and a spline set!\n");return -1;}
    nloops=1;
    nloops1=1;
    if (fdim==3) {
      readitem = strtok(line,"=");
    } else {
      for (i=0; i<fdim-1; i++) {
        if(i==0) {
          readitem = strtok(line,",");
        } else {
          readitem = strtok(NULL,",");
        }
        readitem = strtok(NULL,",");
        strcpy(arSet[i].arIndx,readitem);
        readitem = strtok(NULL,")");
        //strcpy(arSet[i].arSet,readitem);
        //printf("set %s\n",readitem);
        for (i4=0; i4<nset; i4++) if(strcmp(readitem,ha_set[i4].setname)==0) {
            arSet[i].setid=i4;
            //arSet[i].SetSize=ha_set[i4].size;
            //arSet[i].SetBegAdd=ha_set[i4].begadd;
            //arSet[i].subsetid=ha_set[i4].subsetid;
            //if(ha_set[i4].subsetid==1) {
            //  arSet[i].SuperSetSize=ha_set[i4].supersetsize;
            //} else {
            //  arSet[i].SuperSetSize=ha_set[i4].size;
            //}
            break;
          }
        if(i<fdim-3)nloops=nloops*ha_set[arSet[i].setid].size;
        if(i>=fdim-3)nloops1=nloops1*ha_set[arSet[i].setid].size;
      }
      readitem = strtok(NULL,"=");
      dcountdim1[fdim-4]=1;
      for (i=fdim-5; i>-1; i--) {
        dcountdim1[i]=ha_set[arSet[i+1].setid].size*dcountdim1[i+1];
      }
    }
    printf("nloop %ld nloop1 %ld\n",nloops,nloops1);
    strcpy(vname,readitem);
    //readitem = strtok(vname,"(");
        index=ncof-1;
        //printf("vname %s ncof %ld\n",vname,ncof);
        p=strtok(vname,"(");
        //printf("vname %s\n",ha_cof[index].cofname);
        do {
          //printf("ncof %ld indx %ld vname %s\n",ncof,index,ha_cof[index].cofname);
          if (strcmp(ha_cof[index].cofname,p)==0) {
            //printf("vname %s\n",p);
            vvar=index;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
   //printf("read %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){varantidim[l]=0;varsubset[l]=0;varsupsetid[l]=0;}
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  varantidim[dcount]=ha_cof[index].antidims[dcount];
                  vararset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    
    //printf("vname %sloop %ld loop1 %ld ar0 %ld zr1 %ld\n",vname,nloops,nloops1,vararset[0],vararset[1]);
    strcpy(line,linecopy);
    readitem = strtok(line,"=");
    readitem = strtok(NULL,";");
    //printf("read %s\n",readitem);
    i=ha_cgenchf(readitem, ':');
    if(i!=7){printf("Syntax Error in %s!\n",readitem);return -1;}
    strcpy(line,readitem);
    strcat(line,":");
    strcpy(line1,line);
    ///x
    readitem = strtok(line,":");
    strcpy(vname,readitem);
    //readitem = strtok(vname,"(");
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            xvar=index;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
    //printf("read %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){xantidim[l]=0;xsubset[l]=0;xsupsetid[l]=0;}
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  xantidim[dcount]=ha_cof[index].antidims[dcount];
                  xarset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    xsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){xsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
//     if(ha_cof[vvar].size!=ha_cof[xvar].size){
//       printf("Syntax Error!!! %s and %s must have the same size\n",ha_cof[vvar].cofname,ha_cof[xvar].cofname);
//       return 0;
//     }
//     for(l=0;l<ha_cof[vvar].size;l++){
//       if(vararset[l]!=xarset[l]){
//         printf("Syntax Error!!! %s and %s must have the same arguments\n",ha_cof[vvar].cofname,ha_cof[xvar].cofname);
//         return 0;
//       }
//     }
    //printf("vname %sloop %ld loop1 %ld ar0 %ld zr1 %ld\n",vname,nloops,nloops1,xarset[0],xarset[1]);
    ///y
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
    //readitem = strtok(vname,"(");
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            yvar=index;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
    //printf("read %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){yantidim[l]=0;ysubset[l]=0;ysupsetid[l]=0;}
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  yantidim[dcount]=ha_cof[index].antidims[dcount];
                  yarset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    ysubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){ysupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
//     if(ha_cof[vvar].size!=ha_cof[yvar].size){
//       printf("Syntax Error!!! %s and %s must have the same size\n",ha_cof[vvar].cofname,ha_cof[yvar].cofname);
//       return 0;
//     }
//     for(l=0;l<ha_cof[vvar].size;l++){
//       if(vararset[l]!=yarset[l]){
//         printf("Syntax Error!!! %s and %s must have the same arguments\n",ha_cof[vvar].cofname,ha_cof[yvar].cofname);
//         return 0;
//       }
//     }
    
    //printf("vname %sloop %ld loop1 %ld ar0 %ld zr1 %ld\n",vname,nloops,nloops1,yarset[0],yarset[1]);
    ///j
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
    //readitem = strtok(vname,"(");
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            jvar=index;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
    //printf("read %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){jantidim[l]=0;jsubset[l]=0;jsupsetid[l]=0;}
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  jantidim[dcount]=ha_cof[index].antidims[dcount];
                  jarset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    jsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){jsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    if(ha_cof[vvar].size!=ha_cof[jvar].size){
      printf("Syntax Error!!! %s and %s must have the same size\n",ha_cof[vvar].cofname,ha_cof[jvar].cofname);
      return 0;
    }
    for(l=0;l<ha_cof[vvar].size;l++){
      if(vararset[l]!=jarset[l]){
        printf("Syntax Error!!! %s and %s must have the same arguments\n",ha_cof[vvar].cofname,ha_cof[jvar].cofname);
        return 0;
      }
    }
    
    //printf("vname %sloop %ld loop1 %ld ar0 %ld zr1 %ld\n",vname,nloops,nloops1,jarset[0],jarset[1]);

    ///w
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
    //readitem = strtok(vname,"(");
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            wvar=index;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
    //printf("read %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){wantidim[l]=0;wsubset[l]=0;wsupsetid[l]=0;}
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  wantidim[dcount]=ha_cof[index].antidims[dcount];
                  warset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    wsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){wsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    
    //printf("vname %sloop %ld loop1 %ld ar0 %ld zr1 %ld\n",vname,nloops,nloops1,warset[0],warset[1]);

    ///s
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
    //readitem = strtok(vname,"(");
        index=nvar-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_var[index].cofname,p)==0) {
            svar=index;
            if(ha_var[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
        sbegadd=ncofele+ha_var[svar].begadd;
    printf("read %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){santidim[l]=0;ssubset[l]=0;ssupsetid[l]=0;}
            for (dcount=0; dcount<ha_var[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  santidim[dcount]=ha_var[index].antidims[dcount];
                  sarset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_var[index].setid[dcount]].size){
                    ssubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){ssupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    if(ha_cof[vvar].size!=ha_cof[svar].size){
      printf("Syntax Error!!! %s and %s must have the same size\n",ha_cof[vvar].cofname,ha_cof[svar].cofname);
      return 0;
    }
    for(l=0;l<ha_cof[vvar].size;l++){
      if(vararset[l]!=sarset[l]){
        printf("Syntax Error!!! %s and %s must have the same arguments\n",ha_cof[vvar].cofname,ha_cof[svar].cofname);
        return 0;
      }
    }

    ///vlmu
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
    printf("var1 %s\n",vname);
    //readitem = strtok(vname,"(");
        index=ncof-1;
        p=strtok(vname,"(");
    //printf("var2 %s\n",p);
    //printf("var2 %s\n",strtok(NULL,")"));
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            vlmu=index;
            //printf("var2 %s d %d\n",strtok(NULL,")"),ha_var[index].size);
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
    //printf("var %s\n",argu);
            break;
          }
        } while (index--);
        //printf("var3 \n");
        for (l=0; l<MAXVARDIM; l++){mvantidim[l]=0;mvsubset[l]=0;mvsupsetid[l]=0;}
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  mvantidim[dcount]=ha_cof[index].antidims[dcount];
                  mvarset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    mvsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){mvsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
//     if(ha_cof[vvar].size!=ha_cof[xvark].size){
//       printf("Syntax Error!!! %s and %s must have the same size\n",ha_cof[vvar].cofname,ha_cof[xvark].cofname);
//       return 0;
//     }
//     for(l=0;l<ha_cof[vvar].size;l++){
//       if(vararset[l]!=xarset[l]){
//         printf("Syntax Error!!! %s and %s must have the same arguments\n",ha_cof[vvar].cofname,ha_cof[xvark].cofname);
//         return 0;
//       }
//     }
    
    ///dmlu
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
    //readitem = strtok(vname,"(");
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            dlmu=index;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
        for (l=0; l<MAXVARDIM; l++){mdantidim[l]=0;mdsubset[l]=0;mdsupsetid[l]=0;}
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  mdantidim[dcount]=ha_cof[index].antidims[dcount];
                  mdarset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    mdsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){mdsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
        for(l=0;l<4;l++)printf("var2 %lf dmlu %ldbeg %ld\n",ha_cofvar[ha_cof[dlmu].begadd+l].varval,dlmu,ha_cof[dlmu].begadd);
        for(l=10;l<14;l++)printf("var2 %lf\n",ha_cofvar[ha_cof[dlmu].begadd+l].varval);

    ///wlmu
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
    //readitem = strtok(vname,"(");
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            wlmu=index;
            if(ha_var[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
        for (l=0; l<MAXVARDIM; l++){mwantidim[l]=0;mwsubset[l]=0;mwsupsetid[l]=0;}
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  mwantidim[dcount]=ha_cof[index].antidims[dcount];
                  mwarset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    mwsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){mwsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    
    rsize=ha_set[arSet[fdim1-5].setid].size;
    tsize=ha_set[arSet[fdim1-4].setid].size;
    stosize=ha_set[arSet[fdim1-2].setid].size;
    bsize=ha_set[arSet[fdim1-1].setid].size;
    printf("rsize %ld tsize %ld stosize %ld bsize %ld\n",rsize,tsize,stosize,bsize);
    if(bsize!=4)printf("Errors!!! Splinter set with wrong size or in wrong position!!!\n");
    ha_cgetype ***array = (ha_cgetype ***) malloc(sizeof(ha_cgetype **)*rsize);
    for (i = 0; i < rsize; i++)array[i] = (ha_cgetype**) malloc(sizeof(ha_cgetype *)*stosize);
    ha_cgetype vecy2;
    ha_cgetype *matuw= (ha_cgetype *) calloc (rsize*stosize*stosize*bsize,sizeof(ha_cgetype));
    ha_cgetype *matuwold= (ha_cgetype *) calloc (rsize*stosize*stosize*bsize,sizeof(ha_cgetype));
    ha_cgetype *matuv= (ha_cgetype *) calloc (rsize*stosize*stosize,sizeof(ha_cgetype));
    ha_cgetype *matuk= (ha_cgetype *) calloc (rsize*stosize,sizeof(ha_cgetype));
    ha_cgetype *curk= (ha_cgetype *) calloc (rsize,sizeof(ha_cgetype));
    uvadd *curpos= (uvadd *) calloc (rsize,sizeof(uvadd));
    ha_cgetype *vecy= (ha_cgetype *) calloc (rsize*tsize,sizeof(ha_cgetype));
    ha_cgetype *vecx= (ha_cgetype *) calloc (rsize*tsize,sizeof(ha_cgetype));
    ha_cgetype *curx= (ha_cgetype *) calloc (stosize,sizeof(ha_cgetype));
    ha_cgetype *cury= (ha_cgetype *) calloc (stosize,sizeof(ha_cgetype));
    ha_cgetype *curw= (ha_cgetype *) calloc (stosize*bsize,sizeof(ha_cgetype));
    ha_cgetype *vecwf= (ha_cgetype *) calloc (ha_cof[wvar].matsize,sizeof(ha_cgetype));
    ha_cgetype *vecw=NULL,*vecwold=NULL;
    ha_cgetype *vecwoldf= (ha_cgetype *) calloc (ha_cof[wvar].matsize,sizeof(ha_cgetype));
    vecw=vecwf;
    vecwold=vecwoldf;
    //printf("size %ld\n",tsize);
    //strcpy(line1,readitem);
//               if(IsIni&&(IsSplint==1)){
//                 for(i=0;i<nvarele;i++){
//                   if(ha_cgeshock[i].ShockId==1)ha_cgeshock[i].ShockVal=0;//Only mup change in next iter
//                 }
//               }
    if(IsChange==false){
          if((!IsIni)||(IsSplint==2)){
                fout = fopen(filename, "rb");
                if (fout==NULL)printf("Weight file opening error\n");
                freadresult=fread(matuwold,sizeof(ha_cgetype),rsize*stosize*stosize,fout);
                fclose(fout);
          }
          for (l=0; l<rsize; l++){
            arSet[0].indx=l;
            for (i=0; i<stosize; i++){
              arSet[2].indx=i;
              arSet[3].indx=i;//special for uk, only needs 2 set of k
              //dmlu
              l2=0;
              for (dcount=0; dcount<fdim-1; dcount++){
                for(i1=0; i1<ha_cof[dlmu].size; i1++) {
                  if(mvarset[i1]-1==dcount) {
                    if(mvsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[ha_cof[dlmu].setid[i1]].begadd+arSet[dcount].indx].setsh[mdsupsetid[i1]]*mdantidim[i1];
                      //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
                    } else {
                      l2=l2+arSet[dcount].indx*mdantidim[i1];
                    }
                    break;
                  }
                }
              }
              matuk[l*stosize+i]=ha_cofvar[ha_cof[dlmu].begadd+l2].varval;
              printf("k1 %lf v %lf l2 %ld\n",ha_cofvar[ha_cof[dlmu].begadd+l2].varval,ha_cofvar[ha_cof[vlmu].begadd+l2].varval,l2);
              for (j=0; j<stosize; j++){
                arSet[3].indx=j;
                //vlmu
              l2=0;
              for (dcount=0; dcount<fdim-1; dcount++){
                for(i1=0; i1<ha_cof[vlmu].size; i1++) {
                  if(mvarset[i1]-1==dcount) {
                    if(mvsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[ha_cof[vlmu].setid[i1]].begadd+arSet[dcount].indx].setsh[mvsupsetid[i1]]*mvantidim[i1];
                      //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
                    } else {
                      l2=l2+arSet[dcount].indx*mvantidim[i1];
                    }
                    break;
                  }
                }
              }
              matuv[l*stosize*stosize+i*stosize+j]=ha_cofvar[ha_cof[vlmu].begadd+l2].varval;
//               for (i4=0; i4<bsize; i4++){
//                 arSet[4].indx=i4;
//                 //wmlu
//               l2=0;
//               for (dcount=0; dcount<fdim-1; dcount++){
//                 for(i1=0; i1<ha_cof[wlmu].size; i1++) {
//                   if(mwarset[i1]-1==dcount) {
//                     if(mwsubset[i1]==1) {
//                       l2=l2+ha_setele[ha_set[ha_cof[wlmu].setid[i1]].begadd+arSet[dcount].indx].setsh[mwsupsetid[i1]]*mwantidim[i1];
//                       //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
//                     } else {
//                       l2=l2+arSet[dcount].indx*mwantidim[i1];
//                     }
//                     break;
//                   }
//                 }
//               }
//               matuw[l*stosize*stosize*bsize+i*stosize*bsize+j*bsize+i4]=ha_cofvar[ha_cof[vlmu].begadd+l2].varval;
//               }
              }
            }
//           }
//           for (l=0; l<nloops; l++) {
//             i4=l;
//             for (dcount=0; dcount<fdim-3; dcount++) {
//               i3=(uvadd) i4/dcountdim1[dcount];
//               arSet[dcount].indx=i3;
//               i4=i4-i3*dcountdim1[dcount];
//             }
            for(i4=0;i4<tsize;i4++){
              arSet[1].indx=i4;
              //x
              l2=0;
              for (dcount=0; dcount<fdim-1; dcount++){
                for(i1=0; i1<ha_cof[xvar].size; i1++) {
                  if(xarset[i1]-1==dcount) {
                    if(xsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[ha_cof[xvar].setid[i1]].begadd+arSet[dcount].indx].setsh[xsupsetid[i1]]*xantidim[i1];
                      //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
                    } else {
                      l2=l2+arSet[dcount].indx*xantidim[i1];
                    }
                    break;
                  }
                }
              }
              vecx[l*tsize+i4]=ha_cofvar[ha_cof[xvar].begadd+l2].varval;
              //y
              l2=0;
              for (dcount=0; dcount<fdim-2; dcount++){
                for(i1=0; i1<ha_cof[yvar].size; i1++) {
                  if(yarset[i1]-1==dcount) {
                    if(ysubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[ha_cof[yvar].setid[i1]].begadd+arSet[dcount].indx].setsh[ysupsetid[i1]]*yantidim[i1];
                      //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
                    } else {
                      l2=l2+arSet[dcount].indx*yantidim[i1];
                    }
                    break;
                  }
                }
              }
              vecy[l*tsize+i4]=ha_cofvar[ha_cof[yvar].begadd+l2].varval;
              //if(IsSplint==1)vecuy[i]=ha_cofvar[ha_cof[yvar].begadd+l2].varval;
              //else 
//                 vecuy[i4]=ha_cofvar[ha_cof[uyvar].begadd+l2].varval;
//                 if(IsSplint==1)ha_cofvar[ha_cof[uyvar].begadd+l2].varval=vecy[i4];
              printf("var %s l %ld l2 %ld vecy %lf vecx %lf\n",ha_cof[yvar].cofname,l,l2,vecy[l*tsize+i4],vecx[l*tsize+i4]);
            }
          }
            //for (j=0 ; j<(tsize) ; j++)vecy[j]/= vecy[tsize-1];
            //for (j=0 ; j<(tsize) ; j++)vecy[j]/= vecy[tsize-1];
          for (l=0; l<fdim+1; l++)arSet[l].indx=0;
          for (l=0; l<rsize; l++){
            arSet[0].indx=l;
            if(IsSplint==1||IsSplint==2){
              for(i4=0;i4<stosize;i4++){
                l2=l*stosize*stosize+i4*stosize;
                for(i1=0;i1<stosize;i1++)printf("k %lf v %lf\n",matuk[stosize+i1],matuv[l2+i1]);
                spline(matuv+l2,matuk+stosize,0,0,stosize-1,matuw+l*stosize*stosize*bsize+i4*stosize*bsize,laA);
              }
                memcpy(matuwold,matuw,rsize*stosize*stosize*bsize*sizeof(ha_cgetype));
            } else memcpy(matuw,matuwold,rsize*stosize*stosize*bsize*sizeof(ha_cgetype));
            //printf("l %ld uy %s\n",l,ha_cof[uyvar].cofname);
            }
          for (l=0; l<fdim+1; l++)arSet[l].indx=0;
          for (l=0; l<rsize; l++){
            arSet[0].indx=l;
            for(i4=0;i4<tsize;i4++){
              arSet[1].indx=i4;
              l2=0;
              for (dcount=0; dcount<fdim-2; dcount++){
                for(i1=0; i1<ha_cof[jvar].size; i1++) {
                  if(jarset[i1]-1==dcount) {
                    if(jsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[ha_cof[jvar].setid[i1]].begadd+arSet[dcount].indx].setsh[jsupsetid[i1]]*jantidim[i1];
                      //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
                    } else {
                      l2=l2+arSet[dcount].indx*jantidim[i1];
                    }
                    break;
                  }
                }
              }
              temp3=(1+ha_cofvar[ha_cof[jvar].begadd+l2].varval/100)*vecy[l*tsize+i4];
              printf("temp3 %lf\n",temp3);
              for (i1=0; i1<rsize; i1++)curk[i1]=vecy[i1*tsize+i4];
              curk[l]=temp3;
              for (i1=0; i1<rsize; i1++)printf("curk0 %lf curk1 %lf\n",curk[0],vecy[l*tsize+i4]);
              //printf("Shock %lf\n",ha_cofvar[ha_cof[jvar].begadd+l2].varval);
              //New x position
              for (i1=0; i1<rsize; i1++){
              if(matuk[i1*stosize]<matuk[i1*stosize+stosize-1]){
              if(curk[i1]<=matuk[i1*stosize])curpos[i1]=0;//force beg point
              else{
              for(j=1;j<stosize;j++){
                if(j!=stosize-1){
                if(curk[i1]<matuk[i1*stosize+j]&&curk[i1]>=matuk[i1*stosize+j-1]){
                  curpos[i1]=j-1;
                  break;
                }
                }else{
                  curpos[i1]=j-1;
                }
              }
              }
              }else{
              if(curk[i1]>=matuk[i1*stosize])curpos[i1]=0;//force beg point
              else{
              for(j=1;j<stosize;j++){
                if(j!=stosize-1){
                if(curk[i1]>matuk[i1*stosize+j]&&curk[i1]<=matuk[i1*stosize+j-1]){
                  curpos[i1]=j-1;
                  break;
                }
                }else{
                  curpos[i1]=j-1;
                }
              }
              }
              }
              printf("i1 %ld curk %lf cupos %ld\n",i1,curk[i1],curpos[i1]);
              }
              temp3=curk[rsize-1];
              for(j=0;j<stosize;j++){
                l2=l*stosize*stosize*bsize+j*stosize*bsize+curpos[rsize-1];
                curx[j]=matuwold[l2]+matuwold[l2+stosize-1]*temp3+matuwold[l2+2*(stosize-1)]*temp3*temp3+matuwold[l2+3*(stosize-1)]*temp3*temp3*temp3;
              }
              for(j=0;j<stosize;j++)printf("l %ld curx %lf curk %lfcurpos %ld curk0 %lf\n",l,curx[j],matuk[j],curpos[0],curk[0]);
              spline(curx,matuk,0,0,stosize-1,curw,laA);
              temp1=curw[curpos[0]]+curw[curpos[0]+stosize-1]*curk[0]+curw[curpos[0]+2*(stosize-1)]*curk[0]*curk[0]+curw[curpos[0]+3*(stosize-1)]*curk[0]*curk[0]*curk[0];

              l2=0;
              for (dcount=0; dcount<fdim-1; dcount++){
                for(i1=0; i1<ha_cof[vvar].size; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[ha_cof[vvar].setid[i1]].begadd+arSet[dcount].indx].setsh[varsupsetid[i1]]*varantidim[i1];
                      //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
                    } else {
                      l2=l2+arSet[dcount].indx*varantidim[i1];
                    }
                    break;
                  }
                }
              }
                address=ha_cof[vvar].begadd+l2;
                //ha_cofvar[address].varval*=(100+ha_cgeshock[ha_var[svar].begadd+l2].ShockVal)/100;//update
                //ha_cofvar[address].csolpupd=temp1;
//                 ha_cofvar[address].varval=temp1;
                temp2=temp1/ha_cofvar[address].varval;
                temp2=temp2*100-100;
                ha_cgeshock[ha_var[svar].begadd+l2].ShockVal=temp2/subints;//shock, linear only, mup must be exo
                printf("l2 %ld shock %lf var %s temp1 %lf vval %lf shock %lf vvar %lf\n",l2,ha_cgeshock[ha_var[svar].begadd+l2].ShockVal,ha_var[svar].cofname,temp1,ha_cofvar[address].varval,ha_cgeshock[ha_var[svar].begadd+l2].ShockVal,ha_cofvar[ha_cof[vvar].begadd+l2].varval);
              
//               if(IsSplint==1){
//               if(IsIni){
//                 ha_cofvar[ha_cof[vvar].begadd+l2].varval=ha_cofvar[ha_cof[xvar].begadd+l2].varval;
//                 ha_cofvar[ha_cof[vvar].begadd+l2].csolpupd=temp1;
//                 temp2=temp1/ha_cofvar[ha_cof[vvar].begadd+l2].varval;
//                 //printf("temp2 %lf\n",temp2);
//                 temp2=temp2*100-100;
//                 
//                 ha_cgeshock[ha_var[svar].begadd+l2].ShockVal=temp2/subints;//shock, linear only, mup must be exo
//                 //if(temp2<0)temp2=-temp2;
//                 //if(maxerr<temp2)maxerr=temp2;
//                 //printf("sb %ld l2 %ld id %d temp2 %lf temp1 %lf xavr %lf\n",sbegadd,l2,ha_cgeshock[ha_var[svar].begadd+l2].ShockId,ha_cgeshock[ha_var[svar].begadd+l2].ShockVal,temp1,ha_cofvar[ha_cof[vvar].begadd+l2].varval);
//               }else{
//                 //printf("ps %lf varval %lf\n",ha_cgeshock[ha_var[svar].begadd+l2].ShockVal,ha_cofvar[sbegadd+l2].varval);
//                 address=ha_cof[vvar].begadd+l2;
//                 ha_cofvar[address].varval*=(100+ha_cgeshock[ha_var[svar].begadd+l2].ShockVal*subints)/100;//update
//                 ha_cofvar[address].csolpupd=temp1;
//                 temp2=temp1/ha_cofvar[address].varval;
//                 temp2=temp2*100-100;
//                 ha_cgeshock[ha_var[svar].begadd+l2].ShockVal=temp2/subints;//shock, linear only, mup must be exo
//                 //if(temp2<0)temp2=-temp2;
//                 //if(maxerr<temp2)maxerr=temp2;
//                 //printf("sb %ld l2 %ld id %d temp2 %lf temp1 %lf xavr %lf shock %lf\n",sbegadd,l2,ha_cgeshock[ha_var[svar].begadd+l2].ShockId,temp2,temp1,ha_cofvar[address].varval,ha_cofvar[sbegadd+l2].varval);
//               }
//               }else{
//                 if(IsSplint==0){
//                 address=ha_cof[vvar].begadd+l2;
//                 temp2=ha_cofvar[address].varval;
//                 ha_cofvar[address].varval*=(100+ha_cgeshock[ha_var[svar].begadd+l2].ShockVal/nsteps)/100;//update
//                 //ha_cofvar[address].varval=temp1;
//                 //ha_cofvar[address].csolpupd=temp2;
//                 }
//                 if(IsSplint==2){
//                 address=ha_cof[vvar].begadd+l2;
//                 temp2=ha_cofvar[address].varval;
//                 ha_cofvar[address].varval=temp4;//update
//                 //address=ha_cof[xvark].begadd+l2;//must correct before run
//                 ha_cofvar[address].varval=temp5;
//                 }
//               }
            }
//             vecw+=bsize*(tsize-1);
//             vecwold+=bsize*(tsize-1);
          }
//           if(IsSplint==1){
//                 if((fout=fopen(filename, "wb"))==NULL) {
//                   printf("Cannot open file.\n");
//                 }
//                 freadresult=fwrite(matuw, sizeof(ha_cgetype),rsize*stosize*stosize,fout);
//                 fclose(fout);
//           }
    }
    printf("size1 %ld max error in percent %lf\n",tsize,maxerr);
    
    for (i = 0; i < rsize; i++) free(array[i]);
    free(matuk);
    free(matuv);
    free(matuw);
    free(matuwold);
    free(curk);
    free(curpos);
    free(vecy);
    free(vecx);
    free(curx);
    free(cury);
    free(curw);
    free(vecwf);
    free(vecwoldf);
    free(arSet);
  }
  fclose(filehandle);
  return j;
}

uvadd hnew_biupd_bk(PetscInt rank,char *fname,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele,ha_cgeexovar *ha_cgeshock,uvadd nvarele,int laA,uvdim subints,bool IsIni,int IsSplint,int nsteps) {
  FILE * filehandle,*fout;
  char commsyntax[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],linecopy[TABREADLINE];
  char vname[NAMESIZE],sumsyntax[NAMESIZE],argu[NAMESIZE];
  char *readitem=NULL,*p=NULL;
  char filename[1024],j1name[1024];
  uvadd i,i1,i3,i4,l,l2=0,j=0,nsumele,dcountdim1[4*MAXVARDIM],nloops,nloops1,tsize,bsize,vvar,xvar,yvar,jvar,wvar,svar,xvark,uwvar,uyvar,varsize,l2v;//m,
  uvdim fdim,fdim1,dcount,sup;
  bool IsChange=false;
  uvadd index=ncof-1, begadd=0,sbegadd=0,address=0;//,simpl=0;
  bool check10=true;
  ha_cgetype zerodivide=0,temp1,temp2,temp3,maxerr=0,temp4,temp5,sx0,sxn;
  uvadd varantidim[MAXVARDIM],varsubset[MAXVARDIM],vararset[MAXVARDIM],varsupsetid[MAXVARDIM];
  uvadd xantidim[MAXVARDIM],xsubset[MAXVARDIM],xarset[MAXVARDIM],xsupsetid[MAXVARDIM];
  uvadd yantidim[MAXVARDIM],ysubset[MAXVARDIM],yarset[MAXVARDIM],ysupsetid[MAXVARDIM];
  uvadd jantidim[MAXVARDIM],jsubset[MAXVARDIM],jarset[MAXVARDIM],jsupsetid[MAXVARDIM];
  uvadd wantidim[MAXVARDIM],wsubset[MAXVARDIM],warset[MAXVARDIM],wsupsetid[MAXVARDIM];
  uvadd santidim[MAXVARDIM],ssubset[MAXVARDIM],sarset[MAXVARDIM],ssupsetid[MAXVARDIM];
  size_t freadresult;
  //printf("line %s cofvar %f\n",line,ha_cofvar[23].varval);
  ha_cgesetindx *arSet1=NULL;
  hcge_calvars *ha_calvar1= NULL;
  strcpy(commsyntax,"splinter");
      if(rank<10)strcpy(j1name,"000");
      if(rank<100&&rank>9)strcpy(j1name,"00");
      if(rank<1000&&rank>99)strcpy(j1name,"0");
      if(rank>=1000)j1name[0]='\0';
      sprintf(filename, "%d",rank);
      strcat(j1name,filename);
      strcpy(filename,"_biupd");
      strcat(filename,j1name);
      strcat(filename,".bin");
  filehandle = fopen(fname,"r");
  while (ha_cgertabl1(commsyntax,filehandle,line,ha_cofvar,ha_cof,ncof,&zerodivide,TABREADLINE)) {
    IsChange=false;
    if(strstr(line, "(change)")!=NULL) {
      IsChange=true;
      ha_cgefrstr1(line, "(change)", "");
    }
    ha_cgefrstr1(line, commsyntax, "");
    while (ha_cgefrstr(line," ", ""));
    while (ha_cgefrchr(line, '[', '('));
    while (ha_cgefrchr(line, ']', ')'));
    while (ha_cgefrchr(line, '{', '('));
    while (ha_cgefrchr(line, '}', ')'));
    strcpy(linecopy,line);
    //printf("line %s cofvar %f\n",line,ha_cofvar[23].varval);
    //printf("line %s\n",line);
    readitem = strtok(line,"=");//Syntax (linear or change)(other sets):(inter set):(spline set)var_plus()=var():x():x_jump():weight();//
    //printf("read %s\n",readitem);
    fdim=ha_cgenchf(readitem, '(');
    if (fdim==1) {
      fdim=fdim+1;
    }
    fdim1=fdim-1;
    ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (fdim+1,sizeof(ha_cgesetindx));
    //printf("fdim %d\n",fdim);
    if(fdim<3){printf("We should at least have an intertemporal variable and a spline set!\n");return -1;}
    nloops=1;
    nloops1=1;
    if (fdim==3) {
      readitem = strtok(line,"=");
    } else {
      for (i=0; i<fdim-1; i++) {
        if(i==0) {
          readitem = strtok(line,",");
        } else {
          readitem = strtok(NULL,",");
        }
        readitem = strtok(NULL,",");
        strcpy(arSet[i].arIndx,readitem);
        readitem = strtok(NULL,")");
        //strcpy(arSet[i].arSet,readitem);
        //printf("set %s\n",readitem);
        for (i4=0; i4<nset; i4++) if(strcmp(readitem,ha_set[i4].setname)==0) {
            arSet[i].setid=i4;
            //arSet[i].SetSize=ha_set[i4].size;
            //arSet[i].SetBegAdd=ha_set[i4].begadd;
            //arSet[i].subsetid=ha_set[i4].subsetid;
            //if(ha_set[i4].subsetid==1) {
            //  arSet[i].SuperSetSize=ha_set[i4].supersetsize;
            //} else {
            //  arSet[i].SuperSetSize=ha_set[i4].size;
            //}
            break;
          }
        if(i<fdim-3)nloops=nloops*ha_set[arSet[i].setid].size;
        if(i>=fdim-3)nloops1=nloops1*ha_set[arSet[i].setid].size;
      }
      readitem = strtok(NULL,"=");
      dcountdim1[fdim-4]=1;
      for (i=fdim-5; i>-1; i--) {
        dcountdim1[i]=ha_set[arSet[i+1].setid].size*dcountdim1[i+1];
      }
    }
    printf("nloop %ld nloop1 %ld\n",nloops,nloops1);
    strcpy(vname,readitem);
    //readitem = strtok(vname,"(");
        index=ncof-1;
        //printf("vname %s ncof %ld\n",vname,ncof);
        p=strtok(vname,"(");
        //printf("vname %s\n",ha_cof[index].cofname);
        do {
          //printf("ncof %ld indx %ld vname %s\n",ncof,index,ha_cof[index].cofname);
          if (strcmp(ha_cof[index].cofname,p)==0) {
            //printf("vname %s\n",p);
            vvar=index;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
    //printf("read %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){varantidim[l]=0;varsubset[l]=0;varsupsetid[l]=0;}
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  varantidim[dcount]=ha_cof[index].antidims[dcount];
                  vararset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    
    //printf("vname %sloop %ld loop1 %ld ar0 %ld zr1 %ld\n",vname,nloops,nloops1,vararset[0],vararset[1]);
    strcpy(line,linecopy);
    readitem = strtok(line,"=");
    readitem = strtok(NULL,";");
    //printf("read %s\n",readitem);
    i=ha_cgenchf(readitem, ':');
    if(i!=7){printf("Syntax Error in %s!\n",readitem);return -1;}
    strcpy(line,readitem);
    strcat(line,":");
    strcpy(line1,line);
    ///x
    readitem = strtok(line,":");
    strcpy(vname,readitem);
    //readitem = strtok(vname,"(");
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            xvar=index;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
    //printf("read %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){xantidim[l]=0;xsubset[l]=0;xsupsetid[l]=0;}
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  xantidim[dcount]=ha_cof[index].antidims[dcount];
                  xarset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    xsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){xsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    if(ha_cof[vvar].size!=ha_cof[xvar].size){
      printf("Syntax Error!!! %s and %s must have the same size\n",ha_cof[vvar].cofname,ha_cof[xvar].cofname);
      return 0;
    }
    for(l=0;l<ha_cof[vvar].size;l++){
      if(vararset[l]!=xarset[l]){
        printf("Syntax Error!!! %s and %s must have the same arguments\n",ha_cof[vvar].cofname,ha_cof[xvar].cofname);
        return 0;
      }
    }
    //printf("vname %sloop %ld loop1 %ld ar0 %ld zr1 %ld\n",vname,nloops,nloops1,xarset[0],xarset[1]);
    ///y
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
    //readitem = strtok(vname,"(");
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            yvar=index;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
    //printf("read %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){yantidim[l]=0;ysubset[l]=0;ysupsetid[l]=0;}
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  yantidim[dcount]=ha_cof[index].antidims[dcount];
                  yarset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    ysubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){ysupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    if(ha_cof[vvar].size!=ha_cof[yvar].size){
      printf("Syntax Error!!! %s and %s must have the same size\n",ha_cof[vvar].cofname,ha_cof[yvar].cofname);
      return 0;
    }
    for(l=0;l<ha_cof[vvar].size;l++){
      if(vararset[l]!=yarset[l]){
        printf("Syntax Error!!! %s and %s must have the same arguments\n",ha_cof[vvar].cofname,ha_cof[yvar].cofname);
        return 0;
      }
    }
    
    //printf("vname %sloop %ld loop1 %ld ar0 %ld zr1 %ld\n",vname,nloops,nloops1,yarset[0],yarset[1]);
    ///j
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
    //readitem = strtok(vname,"(");
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            jvar=index;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
    //printf("read %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){jantidim[l]=0;jsubset[l]=0;jsupsetid[l]=0;}
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  jantidim[dcount]=ha_cof[index].antidims[dcount];
                  jarset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    jsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){jsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    if(ha_cof[vvar].size!=ha_cof[jvar].size){
      printf("Syntax Error!!! %s and %s must have the same size\n",ha_cof[vvar].cofname,ha_cof[jvar].cofname);
      return 0;
    }
    for(l=0;l<ha_cof[vvar].size;l++){
      if(vararset[l]!=jarset[l]){
        printf("Syntax Error!!! %s and %s must have the same arguments\n",ha_cof[vvar].cofname,ha_cof[jvar].cofname);
        return 0;
      }
    }
    
    //printf("vname %sloop %ld loop1 %ld ar0 %ld zr1 %ld\n",vname,nloops,nloops1,jarset[0],jarset[1]);

    ///w
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
    //readitem = strtok(vname,"(");
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            wvar=index;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
    //printf("read %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){wantidim[l]=0;wsubset[l]=0;wsupsetid[l]=0;}
            for (dcount=0; dcount<ha_cof[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  wantidim[dcount]=ha_cof[index].antidims[dcount];
                  warset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
                    wsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){wsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    
    //printf("vname %sloop %ld loop1 %ld ar0 %ld zr1 %ld\n",vname,nloops,nloops1,warset[0],warset[1]);

    ///s
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
    //readitem = strtok(vname,"(");
        index=nvar-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_var[index].cofname,p)==0) {
            svar=index;
            if(ha_var[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
        sbegadd=ncofele+ha_var[svar].begadd;
    //printf("read %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){santidim[l]=0;ssubset[l]=0;ssupsetid[l]=0;}
            for (dcount=0; dcount<ha_var[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
                  santidim[dcount]=ha_var[index].antidims[dcount];
                  sarset[dcount]=l+1;
                  if (ha_set[arSet[l].setid].size!=ha_set[ha_var[index].setid[dcount]].size){
                    ssubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_var[index].setid[dcount]){ssupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    if(ha_cof[vvar].size!=ha_cof[svar].size){
      printf("Syntax Error!!! %s and %s must have the same size\n",ha_cof[vvar].cofname,ha_cof[svar].cofname);
      return 0;
    }
    for(l=0;l<ha_cof[vvar].size;l++){
      if(vararset[l]!=sarset[l]){
        printf("Syntax Error!!! %s and %s must have the same arguments\n",ha_cof[vvar].cofname,ha_cof[svar].cofname);
        return 0;
      }
    }

    ///xk
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
    printf("var %s\n",vname);
    //readitem = strtok(vname,"(");
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            xvark=index;
            if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
    //printf("read %s\n",argu);
//         for (l=0; l<MAXVARDIM; l++){xantidim[l]=0;xsubset[l]=0;xsupsetid[l]=0;}
//             for (dcount=0; dcount<ha_cof[index].size; dcount++) {
//               if(dcount==0)p=strtok(argu,",");
//               else p=strtok(NULL,",");
//               for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
//                   xantidim[dcount]=ha_cof[index].antidims[dcount];
//                   xarset[dcount]=l+1;
//                   if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
//                     xsubset[dcount]=1;
//                     for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){xsupsetid[dcount]=sup;break;}
//                   }
//                   break;
//                 }
//             }
    if(ha_cof[vvar].size!=ha_cof[xvark].size){
      printf("Syntax Error!!! %s and %s must have the same size\n",ha_cof[vvar].cofname,ha_cof[xvark].cofname);
      return 0;
    }
    for(l=0;l<ha_cof[vvar].size;l++){
      if(vararset[l]!=xarset[l]){
        printf("Syntax Error!!! %s and %s must have the same arguments\n",ha_cof[vvar].cofname,ha_cof[xvark].cofname);
        return 0;
      }
    }
    
    ///uy
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
    //readitem = strtok(vname,"(");
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            uyvar=index;
//             if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
//     //printf("read %s\n",argu);
//         for (l=0; l<MAXVARDIM; l++){yantidim[l]=0;ysubset[l]=0;ysupsetid[l]=0;}
//             for (dcount=0; dcount<ha_cof[index].size; dcount++) {
//               if(dcount==0)p=strtok(argu,",");
//               else p=strtok(NULL,",");
//               for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
//                   yantidim[dcount]=ha_cof[index].antidims[dcount];
//                   yarset[dcount]=l+1;
//                   if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
//                     ysubset[dcount]=1;
//                     for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){ysupsetid[dcount]=sup;break;}
//                   }
//                   break;
//                 }
//             }
//     if(ha_cof[vvar].size!=ha_cof[uyvar].size){
//       printf("Syntax Error!!! %s and %s must have the same size\n",ha_cof[vvar].cofname,ha_cof[uyvar].cofname);
//       return 0;
//     }
//     for(l=0;l<ha_cof[vvar].size;l++){
//       if(vararset[l]!=yarset[l]){
//         printf("Syntax Error!!! %s and %s must have the same arguments\n",ha_cof[vvar].cofname,ha_cof[uyvar].cofname);
//         return 0;
//       }
//     }

    ///uw
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
    //readitem = strtok(vname,"(");
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(ha_cof[index].cofname,p)==0) {
            uwvar=index;
//             if(ha_cof[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
    //printf("read %s\n",argu);
//         for (l=0; l<MAXVARDIM; l++){wantidim[l]=0;wsubset[l]=0;wsupsetid[l]=0;}
//             for (dcount=0; dcount<ha_cof[index].size; dcount++) {
//               if(dcount==0)p=strtok(argu,",");
//               else p=strtok(NULL,",");
//               for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].arIndx,p)==0) {
//                   wantidim[dcount]=ha_cof[index].antidims[dcount];
//                   warset[dcount]=l+1;
//                   if (ha_set[arSet[l].setid].size!=ha_set[ha_cof[index].setid[dcount]].size){
//                     wsubset[dcount]=1;
//                     for(sup=1;sup<MAXSUPSET;sup++)if(ha_set[arSet[l].setid].subsetid[sup]==ha_cof[index].setid[dcount]){wsupsetid[dcount]=sup;break;}
//                   }
//                   break;
//                 }
//             }

    //printf("vname %sloop %ld loop1 %ld ar0 %ld zr1 %ld\n",vname,nloops,nloops1,sarset[0],sarset[1]);
    
    tsize=ha_set[arSet[fdim1-2].setid].size;
    bsize=ha_set[arSet[fdim1-1].setid].size;
    if(bsize!=4)printf("Errors!!! Splinter set with wrong size or in wrong position!!!\n");
    ha_cgetype vecy2;
    ha_cgetype *vecy= (ha_cgetype *) calloc (tsize,sizeof(ha_cgetype));
    ha_cgetype *vecuy= (ha_cgetype *) calloc (tsize,sizeof(ha_cgetype));
(tsize,sizeof(ha_cgetype));
    ha_cgetype *vecx= (ha_cgetype *) calloc (tsize,sizeof(ha_cgetype));
    ha_cgetype *vecwf= (ha_cgetype *) calloc (ha_cof[wvar].matsize,sizeof(ha_cgetype));
    ha_cgetype *vecw=NULL,*vecwold=NULL;
    ha_cgetype *vecwoldf= (ha_cgetype *) calloc (ha_cof[wvar].matsize,sizeof(ha_cgetype));
    vecw=vecwf;
    vecwold=vecwoldf;
    //printf("size %ld\n",tsize);
    //strcpy(line1,readitem);
              if(IsIni&&(IsSplint==1)){
                for(i=0;i<nvarele;i++){
                  if(ha_cgeshock[i].ShockId==1)ha_cgeshock[i].ShockVal=0;//Only mup change in next iter
                }
              }
    if(IsChange==false){
          if((!IsIni)||(IsSplint==2)){
                fout = fopen(filename, "rb");
                if (fout==NULL)printf("Weight file opening error\n");
                freadresult=fread(vecwoldf,sizeof(ha_cgetype),ha_cof[wvar].matsize,fout);
                fclose(fout);
          }
          for (l=0; l<nloops; l++) {
            i4=l;
            for (dcount=0; dcount<fdim-3; dcount++) {
              i3=(uvadd) i4/dcountdim1[dcount];
              arSet[dcount].indx=i3;
              i4=i4-i3*dcountdim1[dcount];
            }
            for(i=0;i<tsize;i++){
              arSet[fdim-3].indx=i;
              //x
              l2=0;
              for (dcount=0; dcount<fdim-2; dcount++){
                for(i1=0; i1<ha_cof[xvar].size; i1++) {
                  if(xarset[i1]-1==dcount) {
                    if(xsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[ha_cof[xvar].setid[i1]].begadd+arSet[dcount].indx].setsh[xsupsetid[i1]]*xantidim[i1];
                      //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
                    } else {
                      l2=l2+arSet[dcount].indx*xantidim[i1];
                    }
                    break;
                  }
                }
              }
              vecx[i]=ha_cofvar[ha_cof[xvar].begadd+l2].varval;
              if(i==0)sx0=ha_cofvar[ha_cof[xvark].begadd+l2].varval;
              if(i==tsize-1)sxn=ha_cofvar[ha_cof[xvark].begadd+l2].varval;
              //y
              l2=0;
              for (dcount=0; dcount<fdim-2; dcount++){
                for(i1=0; i1<ha_cof[yvar].size; i1++) {
                  if(yarset[i1]-1==dcount) {
                    if(ysubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[ha_cof[yvar].setid[i1]].begadd+arSet[dcount].indx].setsh[ysupsetid[i1]]*yantidim[i1];
                      //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
                    } else {
                      l2=l2+arSet[dcount].indx*yantidim[i1];
                    }
                    break;
                  }
                }
              }
              vecy[i]=ha_cofvar[ha_cof[yvar].begadd+l2].varval;
              //if(IsSplint==1)vecuy[i]=ha_cofvar[ha_cof[yvar].begadd+l2].varval;
              //else 
                vecuy[i]=ha_cofvar[ha_cof[uyvar].begadd+l2].varval;
                if(IsSplint==1)ha_cofvar[ha_cof[uyvar].begadd+l2].varval=vecy[i];
              //printf("i %ld l2 %ld vecy %lf\n",i,l2,vecy[i]);
            }
//             for (j=0 ; j<(tsize-1) ; j++){
//               for (i=0 ; i<(tsize-1) ; i++){
//                 if (vecy[i+1] < vecy[i]){
//                   vecy2 = vecy[i];
//                   vecy[i] = vecy[i+1];
//                   vecy[i+1]=vecy2;
//                   vecy2 = vecx[i];
//                   vecx[i] = vecx[i+1];
//                   vecx[i+1]=vecy2;
//                 }
//               }
//             }
//             if(IsSplint==1)for(i=0;i<tsize;i++){
//               arSet[fdim-3].indx=i;
//               //y
//               l2=0;
//               for (dcount=0; dcount<fdim-2; dcount++){
//                 for(i1=0; i1<ha_cof[yvar].size; i1++) {
//                   if(yarset[i1]-1==dcount) {
//                     if(ysubset[i1]==1) {
//                       l2=l2+ha_setele[ha_set[ha_cof[yvar].setid[i1]].begadd+arSet[dcount].indx].setsh[ysupsetid[i1]]*yantidim[i1];
//                       //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
//                     } else {
//                       l2=l2+arSet[dcount].indx*yantidim[i1];
//                     }
//                     break;
//                   }
//                 }
//               }
//               ha_cofvar[ha_cof[uyvar].begadd+l2].varval=vecy[i];
//             }
            for (j=0 ; j<(tsize) ; j++){
              vecy[j]/= vecy[tsize-1];
              vecuy[j]/= vecuy[tsize-1];
            }
            //sx0=0;
            //sxn=0;
            if(IsSplint==1)spline(vecx,vecy,sx0,sxn,tsize-1,vecw,laA);
            else memcpy(vecw,vecwold,bsize*(tsize-1)*sizeof(ha_cgetype));
            //printf("l %ld uy %s\n",l,ha_cof[uyvar].cofname);
            //w
            for(i=0;i<tsize-1;i++){
              arSet[fdim-3].indx=i;//var must be in that orer!!!!
              for(j=0;j<4;j++){
                arSet[fdim-2].indx=j;
              l2=0;
              for (dcount=0; dcount<fdim-1; dcount++){
                for(i1=0; i1<ha_cof[wvar].size; i1++) {
                  if(warset[i1]-1==dcount) {
                    if(wsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[ha_cof[wvar].setid[i1]].begadd+arSet[dcount].indx].setsh[wsupsetid[i1]]*wantidim[i1];
                      //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
                    } else {
                      l2=l2+arSet[dcount].indx*wantidim[i1];
                    }
                    break;
                  }
                }
              }
              //vecwold[j*(tsize-1)+i]=ha_cofvar[ha_cof[wvar].begadd+l2].varval;
              ha_cofvar[ha_cof[wvar].begadd+l2].varval=(ha_floattype)vecw[j*(tsize-1)+i];
              //printf("var %s l2 %ld vi %ld v %lf tem %lf\n",ha_cof[wvar].cofname,l2,j*(tsize-1)+i,vecwold[j*(tsize-1)+i],vecw[j*(tsize-1)+i]);
              }
            }
            for(i=0;i<tsize;i++){
              arSet[fdim-3].indx=i;
              l2=0;
              for (dcount=0; dcount<fdim-2; dcount++){
                for(i1=0; i1<ha_cof[jvar].size; i1++) {
                  if(jarset[i1]-1==dcount) {
                    if(jsubset[i1]==1) {
                      l2=l2+ha_setele[ha_set[ha_cof[jvar].setid[i1]].begadd+arSet[dcount].indx].setsh[jsupsetid[i1]]*jantidim[i1];
                      //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
                    } else {
                      l2=l2+arSet[dcount].indx*jantidim[i1];
                    }
                    break;
                  }
                }
              }
              temp3=(1+ha_cofvar[ha_cof[jvar].begadd+l2].varval/100)*vecy[i];
              //printf("Shock %lf\n",ha_cofvar[ha_cof[jvar].begadd+l2].varval);
              //new x position old weight for error calculation
//               if(!IsIni){
//               temp3=ha_cofvar[ha_cof[vvar].begadd+l2].csolpupd;
//               if(temp3<=vecy[0])j=0;//force beg point
//               else{
//               for(j=1;j<tsize;j++){
//                 if(temp3<=vecy[j]){
//                   j-=1;
//                   break;
//                 }
//               }
//               }
//               if(j>=tsize-1)j=tsize-2;//force end point
//               if(temp3<vecy[0])temp4=vecx[0]+(vecw[j+tsize-1]+2*vecw[j+2*(tsize-1)]*vecy[0]+3*vecw[j+3*(tsize-1)]*vecy[0]*vecy[0])*(temp3-vecy[0]);
//               else if(temp3>vecy[tsize-1])temp4=vecx[tsize-1]+(vecw[j+tsize-1]+2*vecw[j+2*(tsize-1)]*vecy[tsize-1]+3*vecw[j+3*(tsize-1)]*vecy[tsize-1]*vecy[tsize-1])*(temp3-vecy[tsize-1]);
//                 else temp4=vecw[j]+vecw[j+tsize-1]*temp3+vecw[j+2*(tsize-1)]*temp3*temp3+vecw[j+3*(tsize-1)]*temp3*temp3*temp3;
//               //printf("l2 %ld i %ld j %ld x %lf y %lf temp4 %lf j %lf\n",l2,i,j,temp3,vecy[i],temp4,ha_cofvar[ha_cof[jvar].begadd+l2].varval);
//               }
              //New x position
              if(vecuy[0]<=vecuy[tsize-1]){
              if(temp3<=vecuy[0])j=0;//force beg point
              else{
              for(j=1;j<tsize;j++){
                if(j!=tsize-1){
                if(temp3<=vecuy[j]&&temp3>=vecuy[j-1]){
                  j-=1;
                  break;
                }
                if(temp3>=vecuy[j]&&temp3<=vecuy[j-1]){
                  j-=1;
                  break;
                }
                }else{
                if(temp3<=vecuy[j]&&temp3>=vecuy[j-1]){
                  j-=1;
                  break;
                }
                if(temp3>=vecuy[j]&&temp3<=vecuy[j-1]){
                  j-=1;
                  break;
                }
                if(temp3>vecuy[j]){
                  j-=1;
                  break;
                }
                }
              }
              }
              //if(j>=tsize-1)j=tsize-2;//force end point
              }else{
              if(temp3>=vecuy[0])j=0;//force beg point
              else{
              for(j=1;j<tsize;j++){
                if(j!=tsize-1){
                if(temp3<=vecuy[j]&&temp3>=vecuy[j-1]){
                  j-=1;
                  break;
                }
                if(temp3>=vecuy[j]&&temp3<=vecuy[j-1]){
                  j-=1;
                  break;
                }
                }else{
                if(temp3<=vecuy[j]&&temp3>=vecuy[j-1]){
                  j-=1;
                  break;
                }
                if(temp3>=vecuy[j]&&temp3<=vecuy[j-1]){
                  j-=1;
                  break;
                }
                if(temp3<vecuy[j]){
                  j-=1;
                  break;
                }
                }
              }
              }
              }
              //old weight new position
              //if(!IsIni){//redundant lmuw coefficient reading!!!!
                //printf("var %s l2 %ld vi %d v %lf tem %lf\n",ha_cof[wvar].cofname,tsize-1,0,vecwold[tsize-2],vecw[tsize-2]);
              if(j==0){
                temp4=vecx[0]+(vecwold[j+tsize-1]+2*vecwold[j+2*(tsize-1)]*vecuy[0]+3*vecwold[j+3*(tsize-1)]*vecuy[0]*vecuy[0])*(temp3-vecuy[0]);
                //temp4=vecx[0]+((vecx[0]-vecx[1])/(vecuy[0]-vecuy[1]))*(temp3-vecuy[0]);
              } else if(j==tsize-1){
                temp4=vecx[tsize-1]+(vecwold[j+tsize-1]+2*vecwold[j+2*(tsize-1)]*vecuy[tsize-1]+3*vecwold[j+3*(tsize-1)]*vecuy[tsize-1]*vecuy[tsize-1])*(temp3-vecuy[tsize-1]);
                //temp4=vecx[tsize-1]+((vecx[tsize-1]-vecx[tsize-2])/(vecuy[tsize-1]-vecuy[tsize-2]))*(temp3-vecuy[tsize-1]);
              }else {
                  temp4=vecwold[j]+vecwold[j+tsize-1]*temp3+vecwold[j+2*(tsize-1)]*temp3*temp3+vecwold[j+3*(tsize-1)]*temp3*temp3*temp3;
                  //temp4=(temp3-vecy[j])/(vecy[j+1]-vecy[j])*temp4+(vecy[j+1]-temp3)/(vecy[j+1]-vecy[j])*(vecwold[j+1]+vecwold[j+1+tsize-1]*temp3+vecwold[j+1+2*(tsize-1)]*temp3*temp3+vecwold[j+1+3*(tsize-1)]*temp3*temp3*temp3);
                }
              //printf("l2 %ld i %ld j %ld x %lf y %lf temp4 %lf j %lf\n",l2,i,j,temp3,vecy[i],temp4,ha_cofvar[ha_cof[jvar].begadd+l2].varval);
              //}
              //LMUPK
              if(j==0)temp5=vecwold[j+tsize-1]+2*vecwold[j+2*(tsize-1)]*vecuy[0]+3*vecwold[j+3*(tsize-1)]*vecuy[0]*vecuy[0];
              else if(j==tsize-1)temp5=vecwold[j+tsize-1]+2*vecwold[j+2*(tsize-1)]*vecuy[tsize-1]+3*vecwold[j+3*(tsize-1)]*vecuy[tsize-1]*vecuy[tsize-1];
                else {
                  temp5=vecwold[j+tsize-1]+2*vecwold[j+2*(tsize-1)]*temp3+3*vecwold[j+3*(tsize-1)]*temp3*temp3;
                }
              //new position new weight
              if(vecy[0]<=vecy[tsize-1]){
              if(temp3<=vecy[0])j=0;//force beg point
              else{
              for(j=1;j<tsize;j++){
                if(j!=tsize-1){
                if(temp3<=vecy[j]&&temp3>=vecy[j-1]){
                  j-=1;
                  break;
                }
                if(temp3>=vecy[j]&&temp3<=vecy[j-1]){
                  j-=1;
                  break;
                }
                }else{
                if(temp3<=vecy[j]&&temp3>=vecy[j-1]){
                  j-=1;
                  break;
                }
                if(temp3>=vecy[j]&&temp3<=vecy[j-1]){
                  j-=1;
                  break;
                }
                if(temp3>vecy[j]){
                  j-=1;
                  break;
                }
                }
              }
              }
              //if(j>=tsize-1)j=tsize-2;//force end point
              }else{
              if(temp3>=vecy[0])j=0;//force beg point
              else{
              for(j=1;j<tsize;j++){
                if(j!=tsize-1){
                if(temp3<=vecy[j]&&temp3>=vecy[j-1]){
                  j-=1;
                  break;
                }
                if(temp3>=vecy[j]&&temp3<=vecy[j-1]){
                  j-=1;
                  break;
                }
                }else{
                if(temp3<=vecy[j]&&temp3>=vecy[j-1]){
                  j-=1;
                  break;
                }
                if(temp3>=vecy[j]&&temp3<=vecy[j-1]){
                  j-=1;
                  break;
                }
                if(temp3<vecy[j]){
                  j-=1;
                  break;
                }
                }
              }
              }
              }
              if(j==0){
                temp1=vecx[0]+(vecw[j+tsize-1]+2*vecw[j+2*(tsize-1)]*vecy[0]+3*vecw[j+3*(tsize-1)]*vecy[0]*vecy[0])*(temp3-vecy[0]);
                //temp1=vecx[0]+((vecx[0]-vecx[1])/(vecy[0]-vecy[1]))*(temp3-vecy[0]);
              }else if(j==tsize-1){
                temp1=vecx[tsize-1]+(vecw[j+tsize-1]+2*vecw[j+2*(tsize-1)]*vecy[tsize-1]+3*vecw[j+3*(tsize-1)]*vecy[tsize-1]*vecy[tsize-1])*(temp3-vecy[tsize-1]);
                //temp1=vecx[tsize-1]+((vecx[tsize-1]-vecx[tsize-2])/(vecy[tsize-1]-vecy[tsize-2]))*(temp3-vecy[tsize-1]);
              }else {
                  temp1=vecw[j]+vecw[j+tsize-1]*temp3+vecw[j+2*(tsize-1)]*temp3*temp3+vecw[j+3*(tsize-1)]*temp3*temp3*temp3;
                  //temp1=(temp3-vecy[j])/(vecy[j+1]-vecy[j])*temp1+(vecy[j+1]-temp3)/(vecy[j+1]-vecy[j])*(vecw[j+1]+vecw[j+1+tsize-1]*temp3+vecw[j+1+2*(tsize-1)]*temp3*temp3+vecw[j+1+3*(tsize-1)]*temp3*temp3*temp3);
                }
              //printf("i %ld j %ld x %lf y %lf temp1 %lf j %lf\n",i,j,temp3,vecy[i],temp1,ha_cofvar[ha_cof[jvar].begadd+l2].varval);

              if(!IsIni){
                temp2=temp1/temp4;
                temp2=temp2*100-100;
                if(temp2<0)temp2=-temp2;
                if(maxerr<temp2)maxerr=temp2;
                printf("j %ld t1 %lf t2 %lf t3 %lf t4 %lf,w1 %lf w2 %lfw3 %lf w4 %lf\n",j,i*1.0,vecy[i],vecwold[tsize-1]+2*vecwold[2*(tsize-1)]*vecuy[0]+3*vecwold[3*(tsize-1)]*vecuy[0]*vecuy[0],temp5,vecx[i],temp4,temp1,temp2);//if(i>0)(vecw[i+tsize-2]+2*vecw[i+2*(tsize-1)-1]*vecuy[i-1]+3*vecw[i+3*(tsize-1)-1]*vecuy[i-1]*vecuy[i-1]) (vecw[i+tsize-2]+2*vecw[i+2*(tsize-1)-1]*vecuy[i]+3*vecw[i+3*(tsize-1)-1]*vecuy[i]*vecuy[i])
              }
//               l2=0;
//               for (dcount=0; dcount<fdim-2; dcount++){
//                 for(i1=0; i1<ha_cof[vvar].size; i1++) {
//                   if(vararset[i1]-1==dcount) {
//                     if(varsubset[i1]==1) {
//                       l2=l2+ha_setele[ha_set[ha_cof[vvar].setid[i1]].begadd+arSet[dcount].indx].setsh[varsupsetid[i1]]*varantidim[i1];
//                       //printf("supid %d anti %d\n",varsupsetid[i1],varantidim[i1]);
//                     } else {
//                       l2=l2+arSet[dcount].indx*varantidim[i1];
//                     }
//                     break;
//                   }
//                 }
//               }
              if(IsSplint==1){
              if(IsIni){
                ha_cofvar[ha_cof[vvar].begadd+l2].varval=ha_cofvar[ha_cof[xvar].begadd+l2].varval;
                ha_cofvar[ha_cof[vvar].begadd+l2].csolpupd=temp1;
                temp2=temp1/ha_cofvar[ha_cof[vvar].begadd+l2].varval;
                //printf("temp2 %lf\n",temp2);
                temp2=temp2*100-100;
                
                ha_cgeshock[ha_var[svar].begadd+l2].ShockVal=temp2/subints;//shock, linear only, mup must be exo
                //if(temp2<0)temp2=-temp2;
                //if(maxerr<temp2)maxerr=temp2;
                //printf("sb %ld l2 %ld id %d temp2 %lf temp1 %lf xavr %lf\n",sbegadd,l2,ha_cgeshock[ha_var[svar].begadd+l2].ShockId,ha_cgeshock[ha_var[svar].begadd+l2].ShockVal,temp1,ha_cofvar[ha_cof[vvar].begadd+l2].varval);
              }else{
                //printf("ps %lf varval %lf\n",ha_cgeshock[ha_var[svar].begadd+l2].ShockVal,ha_cofvar[sbegadd+l2].varval);
                address=ha_cof[vvar].begadd+l2;
                ha_cofvar[address].varval*=(100+ha_cgeshock[ha_var[svar].begadd+l2].ShockVal*subints)/100;//update
                ha_cofvar[address].csolpupd=temp1;
                temp2=temp1/ha_cofvar[address].varval;
                temp2=temp2*100-100;
                ha_cgeshock[ha_var[svar].begadd+l2].ShockVal=temp2/subints;//shock, linear only, mup must be exo
                //if(temp2<0)temp2=-temp2;
                //if(maxerr<temp2)maxerr=temp2;
                //printf("sb %ld l2 %ld id %d temp2 %lf temp1 %lf xavr %lf shock %lf\n",sbegadd,l2,ha_cgeshock[ha_var[svar].begadd+l2].ShockId,temp2,temp1,ha_cofvar[address].varval,ha_cofvar[sbegadd+l2].varval);
              }
              }else{
                if(IsSplint==0){
                address=ha_cof[vvar].begadd+l2;
                temp2=ha_cofvar[address].varval;
                ha_cofvar[address].varval*=(100+ha_cgeshock[ha_var[svar].begadd+l2].ShockVal/nsteps)/100;//update
                //ha_cofvar[address].varval=temp1;
                //ha_cofvar[address].csolpupd=temp2;
                }
                if(IsSplint==2){
                address=ha_cof[vvar].begadd+l2;
                temp2=ha_cofvar[address].varval;
                ha_cofvar[address].varval=temp4;//update
                address=ha_cof[xvark].begadd+l2;
                ha_cofvar[address].varval=temp5;
                }
              }
            }
            vecw+=bsize*(tsize-1);
            vecwold+=bsize*(tsize-1);
          }
          if(IsSplint==1){
                if((fout=fopen(filename, "wb"))==NULL) {
                  printf("Cannot open file.\n");
                }
                freadresult=fwrite(vecwf, sizeof(ha_cgetype),ha_cof[wvar].matsize,fout);
                fclose(fout);
          }
    }
    printf("size1 %ld max error in percent %lf\n",tsize,maxerr);
    
    free(vecy);
    free(vecuy);
    free(vecx);
    free(vecwf);
    free(vecwoldf);
    free(arSet);
  }
  fclose(filehandle);
  return j;
}
