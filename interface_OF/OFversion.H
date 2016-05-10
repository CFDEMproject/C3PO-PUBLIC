#if(CFDEMWMPROJECTVERSION == 24)
    #define version24
#elif(CFDEMWMPROJECTVERSION == 30)
    #define version30
#endif

// features of 2.4Dev work also in Dev
#if defined(version30)
    #define version24Dev
#endif

// basically use 24x settings + some dev features (e.g. new turbulence model structure)
#if defined(version24Dev)
    #define version24
#endif

// features of 2.4 work also in 2.3
#if defined(version24)
    #define version23
#endif
