#include "array.h"

table* makeAr(int x, int y, int s)
{
     int i=y+1;
     table *t = new table[i];

     *t = (table)i;
     while (--i) t[i]=new char[s*x];

     return t+1;
}

void delAr(table* t)
{
     t--;
     int i=(int)(*t);

     while (--i) delete[]t[i];

     delete[]t;
}
