//#include <ha_cgefparse.h>
//#include <ha_cgetab.h>
#include <ha_cgeglobal.h>
int ha_chrfrall(char *line, int finditem, int replitem)
{
  long int i=0;
  while (line[i]!='\0')
  {
    //printf("i %d\n",i);
    if (line[i]==finditem) {
      //printf("find %d repl %d\n",finditem,replitem);
      line[i]=replitem;
    }
    i++;
  }
  //printf("i %ld\n",i);
  return 0;
}

char* ha_revstrpbrk(char *line, char *finditems)
{
  int count1=0,j=0,l=0,i=0;
  while (line[count1] != '\0') count1++;
  while (finditems[l] != '\0') l++;
  for (i=count1-1; i>-1; i--)
  {
    for(j=0;j<l;j++)if(line[i]==finditems[j]) return &line[i];
  }
  return NULL;
}


int ha_cgerevfind(char *line, char *finditem)
{
  int count1=0,j=0,l=0,i=0;
  char *t1,*t2;
  while (line[count1] != '\0') count1++;
  while (finditem[l] != '\0') l++;
  for (i=count1-1; i>-1; i--)
  {
    t1 = &line[i];
    t2=finditem+l-1;
    j=0;
    while (j<l && tolower((int)*t2)==tolower((int)*t1))
    {
      t1--;
      t2--;
      j++;
    }
    if (j==l) return i;
  }
  return -1;
}

int ha_cgenchf(char *line, int finditem)
{
  int i=0,j=0;
  while (line[i]!='\0')
  {
    if (line[i]==finditem) j++;
    i++;
  }
  return j;
}

int ha_cgenfind(char *line, char *finditem)
{
  int i=0, count1=0,j=0,l=0,q,count2=0;
  char *t1,*t2;
  while (line[count1] != '\0') count1++;
  while (finditem[l] != '\0') l++;
  q=0;
aabb:
  for (i=q; i<count1; i++)
  {
    t1 = &line[i];
    t2=&finditem[0];
    j=0;
    while (j<l && tolower((int)*t2)==tolower((int)*t1))
    {
      t1++;
      t2++;
      j++;
    }
    if (j==l)
    {
      q=i+j;
      count2++;
      goto aabb;
    }
  }
  return count2;
}

char* ha_cgeeqfind(char *line, int finditem)
{
  int count1=0,j=0,i=0,j1,j2;
  while (line[count1] != '\0') count1++;
  for (i=count1-1; i>-1; i--)
  {
    if(line[i]==finditem){
      //printf("i %d\n",i);
      j1=0;j2=0;
      for(j=1;j<count1-i;j++){
        if(line[i+j]==')')j1++;
        if(line[i+j]=='(')j2++;
      }
      //printf("i %d j1 %d j2 %d\n",i,j1,j2);
      if(j1==j2){
        //printf("i %d\n",i);
        return &line[i];
      }
    }
  }
  return NULL;
}
