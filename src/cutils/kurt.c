#include<stdio.h>
#include<stdlib.h>
#include<math.h>


int main(int argc,char *argv[])
{
  int julian,yyyy=1582,mm=10,dd=15;

  julian = julday(mm,dd,yyyy);

  fprintf(stderr,"%d %d %d %d\n",julian,yyyy,mm,dd);

  caldat(julian,&mm,&dd,&yyyy);

  fprintf(stderr,"%d %d %d %d\n",julian,yyyy,mm,dd);


  return 0;
}
