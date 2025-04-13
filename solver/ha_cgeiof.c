//#include <stdio.h>
//#include <string.h>
//#include <stdlib.h>
//#include <ha_cgeiof.h>
#include <ha_cgeglobal.h>

//typedef float ha_cgetype;

int hcge_niodata(char *fname,char *comsyntax) {
  FILE * filehandle;
  char line[TABREADLINE]="\0";
  //char *commsyntax="iodata";
  int j=0;
  filehandle = fopen(fname,"r");
  if(filehandle==NULL){
    printf("Error!!! No such %s file!\n",fname);
    return -1;
  }
  while (ha_cgertabl(comsyntax,filehandle,line,TABREADLINE)) {
    printf("Com %s file %s\n",comsyntax,line);
    j++;
  }
  fclose(filehandle);
  return j;
}

int ha_cgerdvar1(char *varname, char *filename,uvdim *vsize, char *longname,uvdim *d1) {
  FILE * filehandle;
  char line[TABREADLINE+1],linecopy[TABREADLINE+1];
  //char *varname=vname,*filename=fname;
  uvdim nlength=0,nlength1=0,vsizein=0,din1=0;
  int succ=0;
  char *readitem=NULL;
  //printf("varname %s\t%d",varname,nlength);
  while (varname[nlength] != '\0') nlength++;
  //printf("varname %s\t%d\n",varname,nlength);
  filehandle = fopen(filename,"r");
  if(filehandle==NULL){
    printf("Error!!! No such %s file!\n",filename);
    return -1;
  }

  while (fgets(line,TABREADLINE,filehandle)) {
    strcpy(linecopy,line);
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    if (readitem != NULL) {
      nlength1=0;
      while (readitem[nlength1] != '\0'&&readitem[nlength1] != ' ') nlength1++;
      //printf("readname %s\t%d\n",readitem,nlength1);
      if (nlength1==nlength&&strncmp(readitem,varname,nlength) == 0) {
        succ=1;
        readitem = strtok(line," ");
        din1=atoi(readitem);//strtol(readitem,NULL,10);//sscanf(readitem, "%d", &din1);
        readitem = strtok(NULL," ");
        readitem = strtok(NULL," ");
        readitem = strtok(NULL," ");
        vsizein=atoi(readitem);//sscanf(readitem, "%d", &vsizein);
        readitem = strtok(linecopy,"\"");
        readitem = strtok(NULL,"\"");
        readitem = strtok(NULL,"\"");
        readitem = strtok(NULL,"\"");
        strcpy(longname,readitem);
      }
    }
  }
  fclose(filehandle);
  *d1=din1;
  *vsize=vsizein;
  if(succ==0)printf("Error!!! No header \"%s\" in %s\n",varname,filename);
  return succ;
}

int ha_cgermvar1(char *varname, char *filename,uvdim d1, ha_cgemvar1 *record) {
  FILE * dfile;
  char line[DATREADLINE],header[NAMESIZE],varnamecpy[NAMESIZE];
  //char *copyvar;
  uvdim nlength=0,nhead=0,reccount = 0, count1=0,i;
  char *readitem=NULL;
  strcpy(varnamecpy,varname);
  while (ha_cgefrstr(varnamecpy," ", ""));
  while (varnamecpy[nlength] != '\0') nlength++;
  dfile = fopen(filename,"r");

  while (fgets(line,DATREADLINE,dfile)) {
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    if (readitem != NULL) {
      strcpy(header,readitem);
      while (ha_cgefrstr(header," ", ""));
      nhead=0;
      while (header[nhead] != '\0') nhead++;
      if(nhead<nlength)nhead=nlength;
      //printf("var %s head %s nhead %d nlength %d\n",varname,header,nhead,nlength);
      if (strncmp(readitem,varnamecpy,nhead) == 0) {
        //printf("header %s\n",readitem);
        while (fgets(line,DATREADLINE,dfile)) {
          //count1=0;
          //copyvar=line;
          //while (*copyvar==' ' || *copyvar=='\t') copyvar++; // skip whitespaces
          //if (*copyvar=='\r') copyvar++;
          while (ha_cgefrstr(line,"\r", ""));
          while (ha_cgefrstr(line,"  ", " "));
          if (line[0]=='\n') count1=1;
          if (count1!=1) {
            readitem = strtok(line,"\n");
            if (readitem != NULL) {
              record[reccount].dim1=reccount;
              strcpy(record[reccount].ch,readitem);
              reccount++;
              if (reccount>=d1) break;
            }
          }

        }
      }

    }
  }
  fclose(dfile);
  return 0;
}
/*
int ha_cgewmvar1(char *vname, char *longname, char *fname, char *wmode, int rsize, int vsize, int d1, ha_cgemvar1 *record)
{
  FILE * filehandle;
  char line[rsize],longname1[rsize];
  int reccount = 0;
  strcpy(longname1,longname);
  filehandle = fopen(fname,wmode);
  sprintf(line, "%d", d1 );
  fwrite(line, sizeof(char), strlen(line), filehandle);
  fwrite(" Strings Length ", sizeof(char), strlen(" Strings Length "), filehandle);
  sprintf(line, "%d", vsize );
  fwrite(line, sizeof(char), strlen(line), filehandle);
  fwrite(" ", sizeof(char), 1 , filehandle);
  fwrite("\"", sizeof(char), 1 , filehandle);
  fwrite(vname, sizeof(char), strlen(vname), filehandle);
  fwrite("\"", sizeof(char), 1 , filehandle);
  fwrite(" LongName ", sizeof(char), strlen(" LongName ") , filehandle);
  fwrite("\"", sizeof(char), 1 , filehandle);
  fwrite(longname, sizeof(char), strlen(longname), filehandle);
  fwrite("\"", sizeof(char), 1 , filehandle);
  fwrite("\n", sizeof(char), 1 , filehandle);
  while (reccount<d1)
  {
    strcpy(line,record[reccount].ch);
    fwrite(line, sizeof(char), strlen(line), filehandle);
    fwrite("\n", sizeof(char), 1 , filehandle);
    reccount++;
  }
  //fwrite("\n", sizeof(char), 1 , filehandle);
  fclose(filehandle);
  return 0;
}*/

int hcge_rcmd(char *filename, int niodata, hcge_iodata *iodata, char *tabfile, char *closure, char *shock) {
  FILE * filehandle;
  char line[TABREADLINE],*readitem,commsyntax[NAMESIZE];
  int j=0,k;
  strcpy(commsyntax,"iodata");
  filehandle = fopen(filename,"r");
  if(filehandle==NULL){
    printf("Error!!! No such %s file!\n",filename);
    return -1;
  }
  while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
    //printf("line %s\n",line);
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    //printf("read %s\n",readitem);
    k=0;
    while (readitem[k]!= '\0') {
      readitem[k]=tolower((int)readitem[k]);
      k++;
    }
    //printf("read %s\n",readitem);
    strcpy(iodata[j].logname,readitem);
    k=0;
    readitem = strtok(NULL,"\"");
    readitem = strtok(NULL,"\"");
    strcpy(iodata[j].filname,readitem);
    j++;
  }
  fclose(filehandle);

  strcpy(commsyntax,"outdata");
  filehandle = fopen(filename,"r");
  while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
    //printf("line %s\n",line);
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    //printf("read %s\n",readitem);
    k=0;
    while (readitem[k]!= '\0') {
      readitem[k]=tolower((int)readitem[k]);
      k++;
    }
    //printf("read %s\n",readitem);
    strcpy(iodata[j].logname,readitem);
    k=0;
    readitem = strtok(NULL,"\"");
    readitem = strtok(NULL,"\"");
    strcpy(iodata[j].filname,readitem);
    j++;
  }
  fclose(filehandle);

  strcpy(commsyntax,"soldata");
  filehandle = fopen(filename,"r");
  while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
    //printf("line %s\n",line);
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    //printf("read %s\n",readitem);
    k=0;
    while (readitem[k]!= '\0') {
      readitem[k]=tolower((int)readitem[k]);
      k++;
    }
    //printf("read %s\n",readitem);
    strcpy(iodata[j].logname,readitem);
    k=0;
    readitem = strtok(NULL,"\"");
    readitem = strtok(NULL,"\"");
    strcpy(iodata[j].filname,readitem);
    j++;
  }
  fclose(filehandle);

  filehandle = fopen(filename,"r");
  strcpy(commsyntax,"tabfile");
  while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    strcpy(tabfile,readitem);
  }
  fclose(filehandle);
  filehandle = fopen(filename,"r");
  strcpy(commsyntax,"closure");
  while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    strcpy(closure,readitem);
  }
  fclose(filehandle);
  filehandle = fopen(filename,"r");
  strcpy(commsyntax,"shock");
  while (ha_cgertabl(commsyntax,filehandle,line,TABREADLINE)) {
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    strcpy(shock,readitem);
  }
  fclose(filehandle);
  return 1;
}

int ha_csumindx(char *formulain) {
  char line[TABREADLINE],finditem[TABREADLINE],replitem[TABREADLINE],newreplitem[TABREADLINE],temp[TABREADLINE],replitem1[TABREADLINE],newreplitem1[TABREADLINE];
  char*readitem,*line1;
  char syntax[]="sum(";
  int nsum,i,j,l,k,k1,k2;
  nsum=hcge_nsum(formulain,syntax);
  strcpy(line,formulain);
  readitem=line;
  //printf("%d readitem %s\n",nsum,readitem);
  l=0;
  for(i=0;i<nsum-1;i++){
    strcpy(line,formulain);
    k=ha_cgefind(readitem,syntax);
    readitem=readitem+k;
    //printf("readitem1 %s\n",readitem);
    k1=ha_cgefind(readitem,",");
    for(j=0;j<k1+1;j++)finditem[j]=readitem[j];
    finditem[j]='\0';
    readitem=readitem+k1;
    k2=ha_cgefind(readitem,finditem);
    //printf("nsum %d k %d k1 %d k2 %d find %s\n",nsum,k,k1,k2,finditem);
    while(k2!=-1) {
      line1=readitem+k2;
      //printf("readitem2 %s\n",line1);
      sprintf(temp, "%d", l);
      strcat(temp,"?");
      l++;
      for(j=4;j<k1;j++){
        replitem[j-4]=finditem[j];
        newreplitem[j-4]=finditem[j];
      }
      //printf("j1 %d temp1 %s ind %s findiem %s temp %s\n",j,replitem,newreplitem,finditem,temp);
      replitem[j-4]='\0';
      newreplitem[j-4]='\0';
      //printf("j2 %d temp1 %s ind %s findiem %s temp %s\n",j,replitem,newreplitem,finditem,temp);
      strcat(newreplitem,temp);
      ha_cgecutsum(line1);
      strcpy(temp,line1);
      //printf("j %d temp1 %s ind %s findiem %s temp %s\n",j,replitem,newreplitem,finditem,temp);
      //printf("temp %s\n",line1);
      replitem1[0]='(';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,")");
      newreplitem1[0]='(';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,")");
      //printf("repl %s n %s\n",replitem1,newreplitem1);
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]='(';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"]");
      newreplitem1[0]='(';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"]");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]='(';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"}");
      newreplitem1[0]='(';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"}");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]='(';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,",");
      newreplitem1[0]='(';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,",");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]='(';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1," ");
      newreplitem1[0]='(';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1," ");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));


      replitem1[0]='[';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,")");
      newreplitem1[0]='[';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,")");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]='[';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"]");
      newreplitem1[0]='[';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"]");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]='[';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"}");
      newreplitem1[0]='[';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"}");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]='[';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,",");
      newreplitem1[0]='[';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,",");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]='[';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1," ");
      newreplitem1[0]='[';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1," ");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]='{';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,")");
      newreplitem1[0]='{';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,")");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]='{';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"]");
      newreplitem1[0]='{';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"]");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]='{';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"}");
      newreplitem1[0]='{';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"}");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]='{';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,",");
      newreplitem1[0]='{';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,",");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]='{';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1," ");
      newreplitem1[0]='{';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1," ");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]=' ';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,")");
      newreplitem1[0]=' ';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,")");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]=' ';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"]");
      newreplitem1[0]=' ';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"]");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]=' ';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"}");
      newreplitem1[0]=' ';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"}");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]=' ';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,",");
      newreplitem1[0]=' ';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,",");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]=' ';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1," ");
      newreplitem1[0]=' ';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1," ");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]=',';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,")");
      newreplitem1[0]=',';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,")");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]=',';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"]");
      newreplitem1[0]=',';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"]");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]=',';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"}");
      newreplitem1[0]=',';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"}");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]=',';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,",");
      newreplitem1[0]=',';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,",");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      replitem1[0]=',';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1," ");
      newreplitem1[0]=',';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1," ");
      while(ha_cgefrstr(line1,replitem1,newreplitem1));

      //while(ha_cgefrstr(formulain,temp,line1));
      ha_cgefrstr(formulain,temp,line1);
      strcpy(line,formulain);
      //printf("form %s\n",formulain);
      //printf("finfo %s find %s rpl %s\n",formulain,replitem1,newreplitem1);
      //printf("read4 %s find %s\n",readitem,finditem);
      k2=ha_cgefind(readitem,finditem);
    }
  }
  return 0;
}


int hcge_wtab(char *filename, char *newtabfile) {
  //printf("here\n");
  FILE * filehandle,*fout;
  char line[TABREADLINE]="\0",line1[TABREADLINE],line2[TABREADLINE],indx[NAMESIZE],indx1[NAMESIZE],indx2[NAMESIZE],*readitem,*readitem1,commsyntax[NAMESIZE],readline[TABREADLINE],readline1[TABREADLINE],*n,newtabfile1[TABREADLINE];
  char setname[NAMESIZE],newset[NAMESIZE],varname[NAMESIZE],*n1,setelement[TABREADLINE];//,*ne,*np;//,*n2;
  char msetele[MAXVARDIM][NAMESIZE],msetsymb[MAXVARDIM][NAMESIZE],mset[MAXVARDIM][NAMESIZE];
  filehandle = fopen(filename,"r");
  if(filehandle==NULL){
    printf("Error!!! No such %s file!",filename);
    return -1;
  }
  int check,i1,i2,i,setindx,varindx,l1,l2,l3,l4,k1,k2,j,j1,j2;//,necheck,npcheck;//,j;,check1
  strcpy(newtabfile1,newtabfile);
  ha_cgefrstr(newtabfile1,".","1.");
  fout = fopen(newtabfile1,"w");
  readline[0]='\0';
  commsyntax[0]='\0';
  while (fgets(line,TABLINESIZE,filehandle)) {
    //printf("line: %s strlen %d\n",line,strlen(line));
    ha_cgefrstr(line,"\v"," ");
    ha_cgefrstr(line,"![[!","\v");
    ha_cgefrstr(line,"!]]!","\v");
    strcat(readline,line);
    n=strrchr(line,';');
    i1=ha_cgenchf(readline,'!');
    //printf("n %s i1 %d i2 %d\n",readline,i1,i2);
    while (i1>1) {
      while (ha_cgedrcmt(readline,"!"));
      i1-=2;
    }
    //printf("n %s i1 %d i2 %d\n",readline,i1,i2);
    i2=ha_cgenchf(readline,'\v');
    while (i2>1) {
      while (ha_cgedrcmt(readline,"\v"));
      i2-=2;
    }
    //printf("n %s i1 %d i2 %d r len %d\n",readline,i1,i2,strlen(readline));
    if (n!=NULL&&i1==0&&i2==0) {
      check=0;
      //printf("readline: %s strlen %d\n",readline,strlen(readline));
      while (ha_cgedrcmt(readline,"#"));
      while (ha_cgefrstr(readline,"\n", " "));
      while (ha_cgefrstr(readline,"\r", " "));
      while (ha_cgefrstr(readline,"  ", " "));
      while (ha_cgefrstr(readline,", ", ","));
      while (ha_cgefrstr(readline," ,", ","));
      while (ha_cgefrstr(readline,"( ", "("));
      while (ha_cgefrstr(readline," )", ")"));
      k1=0;
      k2=0;
      //printf("read %s\n",readline);
      while (readline[k1]!= '\0') {
        if(readline[k1]=='\"') {
          if(k2==0) k2=1;
          else k2=0;
        } else {
          if(k2==0)readline[k1]=tolower((int)readline[k1]);
        }
        k1++;
      }
      ha_cgefrstr1(readline,"set(", "set (");
      ha_cgefrstr1(readline,"set[", "set [");
      ha_cgefrstr1(readline,"set{", "set {");
      if(ha_cgefind(readline,"set ")==1||ha_cgefind(readline,"set ")==0) {
        strcpy(commsyntax,"set");
        //printf("read0 %s\n",readline);
        ha_cgefrstr(readline,"(", " (");
        while (ha_cgefrstr(readline,"  ", " "));
        //printf("read %s\n",readline);
        //printf("i %d\n",ha_cgefind(readline,"subset "));
        check=1;
      }
      if(check==0){
          while (ha_cgefrstr(readline,"[", "("));
          while (ha_cgefrstr(readline,"]", ")"));
          while (ha_cgefrstr(readline,"{", "("));
          while (ha_cgefrstr(readline,"}", ")"));
      }
      ha_cgefrstr1(readline,"equation(", "equation (");
      ha_cgefrstr1(readline,"formula(", "formula (");
      ha_cgefrstr1(readline,"coefficient(", "coefficient (");
      ha_cgefrstr1(readline,"variable(", "variable (");
      ha_cgefrstr1(readline,"update(", "update (");
      ha_cgefrstr1(readline,"splinter(", "splinter (");
      ha_cgefrstr1(readline,"read(", "read (");
      ha_cgefrstr1(readline,"write(", "write (");
      ha_cgefrstr1(readline,"zerodivide(", "zerodivide (");
      //      if(ha_cgefind(readline,"set ")>-1) {
//        strcpy(commsyntax,"set");
//        check=1;
//      }
      //printf("read %s\n",readline);
      if(ha_cgefind(readline,"subset ")==1||ha_cgefind(readline,"subset ")==0) {
        strcpy(commsyntax,"subset");
        //printf("i %d\n",ha_cgefind(readline,"subset "));
        check=1;
      }
      if(ha_cgefind(readline,"file ")==1||ha_cgefind(readline,"file ")==0) {
        strcpy(commsyntax,"file");
        //printf("i %d\n",ha_cgefind(readline,"file "));
        check=1;
      }
      if (ha_cgefind(readline,"coefficient ")==1||ha_cgefind(readline,"coefficient ")==0) {//if (ha_cgefind(readline,"coefficient ")>-1) {
        strcpy(commsyntax,"coefficient");
        check=1;
      }
      if (ha_cgefind(readline,"variable ")==1||ha_cgefind(readline,"variable ")==0) {//if (ha_cgefind(readline,"variable ")>-1) {
        strcpy(commsyntax,"variable");
        check=1;
      }
      if((ha_cgefind(readline,"read ")==1||ha_cgefind(readline,"read ")==0)&&check==0&&ha_cgefind(readline,"read elements")==-1) {
        strcpy(commsyntax,"read");
        //printf("read %s\n",readline);
        check=1;
      }
      if (ha_cgefind(readline,"formula ")==1||ha_cgefind(readline,"formula ")==0) {
        strcpy(commsyntax,"formula");
        check=1;
      }
      if (ha_cgefind(readline,"equation ")==1||ha_cgefind(readline,"equation ")==0) {//if (ha_cgefind(readline,"equation ")>-1) {
        strcpy(commsyntax,"equation");
        //printf("read %s\n",readline);
        check=1;
      }

      if(ha_cgefind(readline,"update ")==1||ha_cgefind(readline,"update ")==0) {//if(ha_cgefind(readline,"update ")>-1) {
        strcpy(commsyntax,"update");
        check=1;
      }
      if(ha_cgefind(readline,"splinter ")==1||ha_cgefind(readline,"splinter ")==0) {//if(ha_cgefind(readline,"splinter ")>-1) {
        strcpy(commsyntax,"splinter");
        check=1;
      }

      if(ha_cgefind(readline,"zerodivide ")==1||ha_cgefind(readline,"zerodivide ")==0) {//if(ha_cgefind(readline,"zerodivide ")>-1) {
        strcpy(commsyntax,"zerodivide");
        check=1;
      }

      if(ha_cgefind(readline,"write ")==1||ha_cgefind(readline,"write ")==0) {//if(ha_cgefind(readline,"write ")>-1) {
        strcpy(commsyntax,"write");
        check=1;
      }

      if(ha_cgefind(readline,"assertion ")==1||ha_cgefind(readline,"assertion ")==0) {//if(ha_cgefind(readline,"assertion ")>-1) {
        strcpy(commsyntax,"assertion");
        check=1;
      }
      //printf("com %s\n",commsyntax);
      if (strchr(readline,';')!=NULL) {
        if (check==1) {
          if (readline[0]==' ') fprintf(fout,"%s\n",readline+1);
          else fprintf(fout,"%s\n",readline);
          readline[0]='\0';
        } else {
          if (readline[0]==' ') fprintf(fout,"%s%s\n",commsyntax,readline);
          else fprintf(fout,"%s %s\n",commsyntax,readline);
          readline[0]='\0';
        }
      }
    }
  }
  fclose(filehandle);
  fclose(fout);
  filehandle = fopen(newtabfile1,"r");
  fout = fopen(newtabfile,"w");
  i=0;
  while (fgets(line,TABREADLINE,filehandle)) {
    //printf("line %s\n",line);
    ha_csumindx(line);
    //printf("newline %s\n",line);
    l1=ha_cgefind(line,"equation");
    l2=ha_cgefind(line,"formula");
    l3=ha_cgefind(line,"read");
    l4=ha_cgefind(line,"update");
    if (l3==0&&ha_cgenfind(line,"\"")<3) l3=-1;
    if (l1==0||l2==0||l3==0||l4==0) {
      if (l1==0||l4==0) {
        n=strstr(line,"c_");
        while (n!=NULL) {
          if (line[n-line-1]==' '||line[n-line-1]=='+'||line[n-line-1]=='-'||line[n-line-1]=='*'||line[n-line-1]=='/'||line[n-line-1]=='^'||line[n-line-1]==','||line[n-line-1]=='=') {
            line[n-line]='p';
            n=strstr(line,"c_");
          } else n=strstr(n+1,"c_");
        }
      }
      strcpy(line1,line);
      //printf("line %s\n",line);
      n=strchr(line1,'\"');
      if (n!=NULL) {
        readline[0]='\0';
        readline1[0]='\0';
        //j=0;
        while (n!=NULL) {
          strncpy(line2,line1,n-line1);
          line2[n-line1]='\0';
          n1=strrchr(line2,'(');
          setindx=ha_cgenfind(n1,",");
          varindx=n1-line2+1;
          line2[varindx]='\0';
          //printf("line2 %s\n",line2);
          n1=strrchr(line2,' ');
          if (n1!=NULL) varindx=n1-line2;
          n1=strrchr(line2,',');
          if (n1!=NULL&&(n1-line2)>varindx) varindx=n1-line2;
          n1=strrchr(line2,'*');
          if (n1!=NULL&&(n1-line2)>varindx) varindx=n1-line2;
          n1=strrchr(line2,'/');
          //printf("n1 %d\n",n1-line2);
          if (n1!=NULL&&(n1-line2)>varindx) varindx=n1-line2;
          n1=strrchr(line2,'+');
          if (n1!=NULL&&(n1-line2)>varindx) varindx=n1-line2;
          n1=strrchr(line2,'-');
          if (n1!=NULL&&(n1-line2)>varindx) varindx=n1-line2;
          n1=strrchr(line2,'^');
          if (n1!=NULL&&(n1-line2)>varindx) varindx=n1-line2;
          n1=strrchr(line2,'=');
          if (n1!=NULL&&(n1-line2)>varindx) varindx=n1-line2;
          strcpy(varname,&line2[varindx+1]);
          n1=strrchr(varname,'(');
          varname[n1-varname+1]='\0';
           j1=0;
           j2=0;
           while (varname[j1]!= '\0'){
             if(varname[j1]!= '(')j2++;
             j1++;
           }
           n1=varname;
           if(j2>1){
             j2=0;
             while (j1>-1){
             if(varname[j1]== '(')j2++;
             if(j2==2){
               n1=varname+j1+1;
               break;
             }
             j1--;
             }
           }
          //printf("var %s i %d\n",n1,i);
          hcge_rsetname(newtabfile1,n1,setindx,setname);
          //printf("setname %s\n",setname);
          strcpy(indx,"i");
          sprintf(indx1, "%d",i);
          strcat(indx,indx1);
//           while (strstr(line1,indx)!=NULL) {
//             i++;
//             strcpy(indx,"i");
//             sprintf(indx1, "%d",i);
//             strcat(indx,indx1);
//           }
          strcpy(newset,"sub_");
          strcat(newset,indx);
          //printf("sub set %s\n",newset);
          i1=n-line1;
          readitem=strtok(n,"\"");
          strcpy(setelement,readitem);
          //printf("read %s\n",readitem);
          k1=0;
          while(setelement[k1]!='\0') {
            setelement[k1]=tolower((int)setelement[k1]);
            k1++;
          }
          //printf("linerr %s\n",line);
          //printf("read %s\n",setelement);
          strcpy(indx1,"\"");
          strcat(indx1,readitem);
          strcat(indx1,"\"");
//           for(j1=0;j1<j;j1++)if(strcmp(setelement,msetele[j1])==0&&strcmp(setname,mset[j1])==0){//&&i>78&&i<83){
//             ha_cgefrstr(line,indx1,msetsymb[j1]);
//             break;
//           }
          //if(j1==j){
          strcat(readline,"set ");
          strcat(readline,newset);
          strcat(readline," (");
          strcat(readline,setelement);
          //strcpy(msetele[j],setelement);
          strcat(readline,");\n");
          strcat(readline1,"subset sub_");
          strcat(readline1,indx);
          //strcpy(msetsymb[j],indx);
          strcat(readline1," is subset of ");
          strcat(readline1,setname);
          //strcpy(mset[j],setname);
          strcat(readline1," ;\n");
          ha_cgefrstr(line,indx1,indx);
          //if (l3==0) ha_cgefrstr(line,indx1,indx);
          //else while(ha_cgefrstr(line,indx1,indx));
          //printf("linerr %s\n",line);
          //}
          readitem1=strtok(line," ");
          //printf("read122 %s\n",readitem1);
          strcpy(line2,readitem1);
          strcat(line2," ");
          readitem1=strtok(NULL," ");
          strcpy(indx2,indx);
          strcat(indx2,",");
          //printf("indx %s readitem1 %s\n",indx,readitem1);
          //printf("line2 %s\n",line2);
          if(strstr(readitem1,indx2)==NULL) {
            //printf("line2 %s setname %s\n",line2,indx);
            if(l1==0){//if(l1==0&&strchr(line2,'(')!=NULL) {
              strcat(line2,readitem1);
              strcat(line2," ");
              readitem1=strtok(NULL," ");
            }
            //printf("line2 %s\n",line2);
            //if(j1==j){
            strcat(line2," (all,");
            strcat(line2,indx);
            strcat(line2,",");
            strcat(line2,newset);
            strcat(line2,")");
            //}
            strcat(line2,readitem1);
            strcat(line2," ");
          } else {
            //if(j1==j){
            strcat(line2," (all,");
            strcat(line2,indx);
            strcat(line2,",");
            strcat(line2,newset);
            strcat(line2,")");
            //}
            strcat(line2," ");
            strcat(line2,readitem1);
            //printf("l1 %d line2 %s setname %s\n",l1,line2,indx);
//             if(l1==0&&strchr(line2,'(')!=NULL) {//do nothing
//               readitem1=strtok(NULL," ");
//               strcat(line2," ");
//               strcat(line2,readitem1);
//             }
            //printf("line2 %s read 1 %s\n",line2,readitem1);
          }
          readitem1=strtok(NULL,"\n");
          //printf("line2 %s\n",line2);
          //printf("read1 %s\n",readitem1);
          if(readitem1!=NULL)strcat(line2,readitem1);
          //printf("line2 %s\n",line2);
          strcat(line2,"\n");
          //printf("line2 %s\n",line2);
          strcpy(line,line2);
          //}
          strcpy(line1,line);
          n=strchr(line1,'\"');
          if (l3==0&&ha_cgenfind(line1,"\"")<3) n=NULL;
          i++;
          //j++;
        }
        //printf("read %s read1 %s line2 %s\n",readline,readline1,line2);
        fprintf(fout,"%s%s%s",readline,readline1,line2);
      } else fprintf(fout,"%s",line);
    } else fprintf(fout,"%s",line);
  }
  fclose(filehandle);
  fclose(fout);
  remove(newtabfile1);
  return 1;
}

int hcge_wdata(char *filename, char *newdatlogname, char *newdatfile,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele,hcge_cof *ha_cof,uvadd ncof,uvadd ncofele,hcge_cof *ha_var,uvadd nvar,uvadd nvarele, ha_cgevar *ha_cofvar) {
  //printf("here %s dat %s\n",filename,newdatfile);
  FILE * filehandle,*fout;
  char line[TABREADLINE]="\0",*readline,comsyntax[TABREADLINE],longname[TABREADLINE],datline[DATREADLINE],varname[NAMESIZE],*vname1,header[NAMESIZE],setsize[DATREADLINE],tempname[NAMESIZE];
  filehandle = fopen(filename,"r");
  long int i,n,j,j1,j2,innerloop,outerloop,l,indx;
  long int setindx[MAXVARDIM],antidim[MAXVARDIM];
  strcpy(comsyntax,"to file ");
  strcat(comsyntax,newdatlogname);
  strcat(comsyntax," header \"");
  //printf("comsyntax %s\n",comsyntax);
  n=strlen(comsyntax);
  fout = fopen(newdatfile,"w");
  i=0;
  while (fgets(line,TABREADLINE,filehandle)) {
    if(strncmp(line,"write",5)==0&&strstr(line,comsyntax)!=NULL) {
      if(strstr(line,"(set)")!=NULL) {
        i=ha_cgefind(line," ");
        readline=line+i+1;
        i=ha_cgefind(readline," ");
        readline=readline+i+1;
        i=ha_cgefind(readline," ");
        strncpy(varname,readline,i);
        varname[i]='\0';
        //printf("var %s\n",varname);
        i=ha_cgefind(line,comsyntax);
        readline=line+i+n;
        i=ha_cgefind(readline,"\"");
        strncpy(header,readline,i);
        header[i]='\0';

        i=ha_cgefind(line,"longname \"");
        readline=line+i+10;
        //printf("var %s\n",readline);
        i=ha_cgefind(readline,"\"");
        strncpy(longname,readline,i);
        longname[i]='\0';
        setsize[0]='\0';
        for (i=0; i<nset; i++) {
          if (strcmp(ha_set[i].setname,varname)==0) {
            sprintf(setsize, "%d", ha_set[i].size);
            strcat(setsize," Strings Length 12 Header \"");
            strcat(setsize,header);
            strcat(setsize,"\" LongName \"");
            strcat(setsize,longname);
            strcat(setsize,"\";\n");
            fprintf(fout,"%s",setsize);
            for (j=0; j<ha_set[i].size; j++)fprintf(fout,"%s\n",ha_setele[ha_set[i].begadd+j].setele);
            fprintf(fout,"\n");
          }
        }
      } else {
        //printf("read %s\n",line);
        i=ha_cgefind(line," ");
        readline=line+i+1;
        i=ha_cgefind(readline," ");
        strncpy(varname,readline,i);
        varname[i]='\0';
        //printf("var %s\n",varname);
        i=ha_cgefind(line,comsyntax);
        if(i==-1)continue;
        readline=line+i+n;
        i=ha_cgefind(readline,"\"");
        strncpy(header,readline,i);
        header[i]='\0';

        i=ha_cgefind(line,"longname \"");
        readline=line+i+10;
        //printf("var %s\n",readline);
        i=ha_cgefind(readline,"\"");
        strncpy(longname,readline,i);
        longname[i]='\0';
        //printf("var %s\n",longname);
        //fprintf(fout,"Header \"%s\" LongName \"%s\"\n",header,longname);
        //printf("Header \"%s\" LongName \"%s\"\n",header,longname);
        setsize[0]='\0';
        for (i=0; i<ncof; i++) {
          vname1= strtok(ha_cof[i].cofname,"(");
          //printf("i %d var %s var1 %s\n",i,varname,ha_cof[i].cofname);

          if (strcmp(vname1,varname)==0) {
            if(ha_cof[i].size==0) {
              innerloop=1;
              outerloop=0;
            }
            if(ha_cof[i].size==1) {
              innerloop=ha_set[ha_cof[i].setid[0]].size;
              outerloop=0;
            }
            if(ha_cof[i].size==2) {
              innerloop=ha_set[ha_cof[i].setid[0]].size*ha_set[ha_cof[i].setid[1]].size;
              outerloop=1;
            }
            if(ha_cof[i].size>2) {
              innerloop=ha_set[ha_cof[i].setid[0]].size*ha_set[ha_cof[i].setid[1]].size;
              outerloop=1;
              for(j=2; j<ha_cof[i].size; j++)outerloop*=ha_set[ha_cof[i].setid[j]].size;
            }
            //printf("var %s in %d out %d\n",varname,innerloop,outerloop);
            for(j=0; j<ha_cof[i].size; j++) {
              sprintf(tempname, "%d", ha_set[ha_cof[i].setid[j]].size);
              strcat(setsize,tempname);
              strcat(setsize," ");
              //printf("set %d %s\n",ha_cof[i].antidims[j],ha_set[ha_cof[i].setid[j]].setname);
            }
            if(ha_cof[i].size==0){
              sprintf(tempname, "%d", 1);
              strcat(setsize,tempname);
              strcat(setsize," ");
            }
            strcat(setsize,"Real SpreadSheet Header \"");
            strcat(setsize,header);
            strcat(setsize,"\" LongName \"");
            strcat(setsize,longname);
            strcat(setsize,"\";\n");
            fprintf(fout,"%s",setsize);
            indx=0;
            //printf("var %s\n",setsize);
            //printf("inloop %d s %d outloop %d\n",innerloop,ha_cof[i].size,outerloop);
            //for(j=0; j<ha_cof[i].size; j++)printf("i %d s %d\n",j,ha_set[ha_cof[i].setid[j]].size);
            if(ha_cof[i].size<2) {
              for(j=0; j<innerloop; j++) {
                fprintf(fout,"%f\n",ha_cofvar[ha_cof[i].begadd+j].varval);
              }
              fprintf(fout,"\n");
            } else {
              antidim[2]=1;
              for (l=3; l<ha_cof[i].size; l++){
                antidim[l]=antidim[l-1]*ha_set[ha_cof[i].setid[l-1]].size;
                //printf("var %s l %d anti %d\n",ha_cof[i].cofname,l,antidim[l]);
              }
              for(j1=0; j1<outerloop; j1++) {
                indx=j1;
                for (l=ha_cof[i].size-1; l>1; l--) {
                  setindx[l]=indx/antidim[l];
                  indx-=setindx[l]*antidim[l];
                  //printf("var %s l %d set %d\n",ha_cof[i].cofname,l,setindx[l]);
                }
                for(j=0; j<innerloop; j++) {
                  setindx[0]=j/ha_set[ha_cof[i].setid[1]].size;
                  setindx[1]=j-ha_set[ha_cof[i].setid[1]].size*setindx[0];
                  indx=0;
                  for (l=0; l<ha_cof[i].size; l++)indx+=ha_cof[i].antidims[l]*setindx[l];
                  if(setindx[1]==ha_set[ha_cof[i].setid[1]].size-1){
                    fprintf(fout,"%f\n",ha_cofvar[ha_cof[i].begadd+indx].varval);
                  }else{
                    fprintf(fout,"%f,",ha_cofvar[ha_cof[i].begadd+indx].varval);
                  }
                  //if(strcmp(header,"1BAS")==0)printf("i0 %d i1 %d i2 %d indx %d val \%lf\n",setindx[0],setindx[1],setindx[2],indx,ha_cofvar[ha_cof[i].begadd+indx].varval);
                }
                fprintf(fout,"\n");
              }
            }
            break;
          }
        }
      }
    }
  }
  fclose(filehandle);
  fclose(fout);
  return 1;
}

int hcge_wvar(char *filename, char *newtabfile,hcge_cof *ha_var,uvadd nvar) {
  FILE * filehandle,*fout;
  char line[TABREADLINE+1]="\0",*p;//,nvarname[nvar][NAMESIZE+2],*p;//,line1[DATREADLINE];//,*ne,*np;//,*n2;
  filehandle = fopen(filename,"r");
  uvadd i,n,j,l,l1,linelght;
  int lvar;
  //char *nvarname= (char *) calloc (nvar,sizeof(varname));
  fout = fopen(newtabfile,"w");
  //for (i=0;i<nvar;i++) {
  //if(strncmp(ha_var[i].cofname,"p_",2)==0){
  //strcpy(nvarname[i],ha_var[i].cofname);
  //}
  //else {
  //strcpy(nvarname[i],"p_");
  //strcat(nvarname[i],ha_var[i].cofname);
  //}
  //}
  while (fgets(line,TABREADLINE,filehandle)) {
    if(ha_cgefind(line,"equation ")>-1||ha_cgefind(line,"update ")>-1) {
      //printf("line %s\n",line);
      linelght=strlen(line);
      for (i=0; i<nvar; i++) {
        p=strchr(line,';');
        line[p-line+1]='\n';
        //strcat(line,"\n");
        line[p-line+2]='\0';
        n=ha_cgenfind(line,ha_var[i].cofname);
        //printf("n %d var %s\n",n,ha_var[i].cofname);
        lvar=strlen(ha_var[i].cofname);
        //if(strcmp(ha_var[i].cofname,"x5")==0&&n>0) printf("line %s\n n %d\n",line,n);
        l=0;
        for (j=0; j<n; j++) {
          l1=ha_cgefind(&line[l],ha_var[i].cofname);
          l=l+l1;
          //linelght=strlen(line);
          if(strncmp(ha_var[i].cofname,"p_",2)!=0&&ha_var[i].level_par==false) if(line[l+lvar]==' '||line[l+lvar]=='('||line[l+lvar]=='+'||line[l+lvar]=='-'||line[l+lvar]=='*'||line[l+lvar]=='/'||line[l+lvar]=='^'||line[l+lvar]==']'||line[l+lvar]==','||line[l+lvar]==';'||line[l+lvar]=='=')if(line[l-1]==' '||line[l-1]=='+'||line[l-1]=='-'||line[l-1]=='*'||line[l-1]=='/'||line[l-1]=='^'||line[l-1]=='['||line[l-1]=='('||line[l-1]==','||line[l-1]=='=') {
                memmove(&line[l+2],&line[l],linelght-l);
                line[l]='p';
                line[l+1]='_';
                l=l+2;
                linelght+=2;
              }
          //if(strncmp(ha_var[i].cofname,"c_",2)!=0&&ha_var[i].pchange_real==true)if(line[l-1]==' '||line[l-1]=='+'||line[l-1]=='-'||line[l-1]=='*'||line[l-1]=='/'||line[l-1]=='^'||line[l-1]=='['||line[l-1]=='('||line[l-1]==','){memmove(&line[l+2],&line[l],linelght-l);line[l]='c';line[l+1]='_';l=l+2;}
          l=l+strlen(ha_var[i].cofname);
        }
        //ha_cgefrstr(line,ha_var[i].cofname,nvarname[i]);
      }

      //line[0]='\0';
      //printf("line1 %s\n",line);
    }
    //p=strchr(line,';');
    //line[p-line+1]='\n';
    //line[p-line+2]='\0';
    fprintf(fout,"%s",line);
  }
  fclose(filehandle);
  fclose(fout);
  return 1;
}


int hcge_rsetname(char *filename, char *varname, int indx, char *setname) {
  FILE * filehandle;
  int n,i;
  uvadd lsize;
  lsize=TABREADLINE+1;
  char commsyntax[NAMESIZE],varname1[NAMESIZE+2],indxname[NAMESIZE],line[TABREADLINE+1],line1[TABREADLINE+1],*p,tmp[TABREADLINE+1];//,varname2[NAMESIZE+2],varname3[NAMESIZE+2]
  //printf("var %s\n",varname);
  strcpy(varname1,")");
  strcat(varname1,varname);
  filehandle = fopen(filename,"r");
  strcpy(commsyntax,"variable");
  //printf("var %s\n",varname1);
  while (ha_cgertabl(commsyntax,filehandle,line,lsize)) {
    while (ha_cgefrstr(line," ", ""));
    strcpy(line1,line);
    //printf("line %s\n",line);
    n=ha_cgefind(line,varname1);
    if (n==-1&&((varname[0]=='p'&&varname[1]=='_')||(varname[0]=='c'&&varname[1]=='_'))) {
      strcpy(varname1,")");
      strcat(varname1,varname+2);
      n=ha_cgefind(line,varname1);
    }
    //if(n>0) printf("line %s %d\n",line,indx);
    if (n>-1) {
      p=strtok(line+n,"(");
      //printf("line %s\n",p);
      p=strtok(NULL,")");
      strcpy(tmp,p);
      strcat(tmp,",");
      for (i=0; i<indx+1; i++) {
        if(i==0) p=strtok(tmp,",");
        else p=strtok(NULL,",");
      }
      //if(indx==0) {
      //  p=strtok(tmp,",");
      //}
      //else {
      //  for (i=0; i<indx; i++) {
      //    if(i==0) p=strtok(tmp,",");
      //    else p=strtok(NULL,",");
      //  }
      //  p=strtok(NULL,",");
      //}
      //if(n>0) printf("line %s\n",p);
      strcpy(indxname,",");
      strcat(indxname,p);
      strcpy(line,line1);
      n=ha_cgefind(line,indxname);
      p=strtok(line+n,",");
      p=strtok(NULL,")");
      strcpy(setname,p);
      return 1;
    }
  }
  fclose(filehandle);
  filehandle = fopen(filename,"r");
  strcpy(commsyntax,"coefficient");
  while (ha_cgertabl(commsyntax,filehandle,line,lsize)) {
    while (ha_cgefrstr(line," ", ""));
    strcpy(line1,line);
    //printf("line %s indx %d\n",line,indx);
    n=ha_cgefind(line,varname1);
    if (n==-1&&((varname[0]=='p'&&varname[1]=='_')||(varname[0]=='c'&&varname[1]=='_'))) {
      strcpy(varname1,")");
      strcat(varname1,varname+2);
      n=ha_cgefind(line,varname1);
    }
    //if(n>-1)printf("line %s %d\n",line,n);
    if (n>-1) {
      p=strtok(line+n,"(");
      //printf("line %s\n",p);
      p=strtok(NULL,")");
      strcpy(tmp,p);
      strcat(tmp,",");
      for (i=0; i<indx+1; i++) {
        if(i==0) p=strtok(tmp,",");
        else p=strtok(NULL,",");
      }
      strcpy(indxname,",");
      strcat(indxname,p);
      strcat(indxname,",");
      strcpy(line,line1);
      //printf("line %s\n",indxname);
      n=ha_cgefind(line,indxname);
      p=strtok(line+n,",");
      p=strtok(NULL,")");
      strcpy(setname,p);
      return 1;
    }
  }
  fclose(filehandle);
  return -1;
}

bool intreadCSV(char *fileName,int* vec, int vecCol) {
  FILE * filehandle;
  char line[DATREADLINE],*p;
  int i=0,j;
  filehandle = fopen(fileName,"r");
  j=0;
  while (fgets(line,DATREADLINE,filehandle)) {
    if(j==0){j++;continue;}
    strcat(line,",");
      for(i=0; i<vecCol+1; i++) {
        if(i==0) {
          p = strtok(line,",");
        } else {
          p = strtok(NULL,",");
          vec[i-1]=atoi(p);
        }
      }
    j++;
  }
  fclose(filehandle);
  return true;
}
