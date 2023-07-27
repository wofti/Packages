/* WT: if we set 
#define STATIC 
only we can make many places more thread safe, while
#define STATIC static
should result in the old behavior */
#define STATIC 
