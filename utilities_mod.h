#pragma once
#include <math.h>
#include <string.h>
#include <iostream>
char *itos(int ndigit, int n)
{

   char str[10];
   int i, nmax, nd;
   auto a = std::to_string(32) + "";
   nmax = pow(10, ndigit) - 1;
   // Count no. of zeroes to pad
   if (n <= 9)
   {
      std::string temp = " " + std::to_string(n);
      strcpy(str, temp.c_str());
      nd = ndigit - 1;
   }
   else if (n <= 99)
   {
      std::string temp = "  " + std::to_string(n);
      strcpy(str, temp.c_str());
      nd = ndigit - 2;
   }
   else if (n <= 999)
   {
      std::string temp = "   " + std::to_string(n);
      strcpy(str, temp.c_str());
      nd = ndigit - 3;
   }
   else if (n <= 9999)
   {
      std::string temp = "    " + std::to_string(n);
      strcpy(str, temp.c_str());
      nd = ndigit - 4;
   }
   else
   {
      nd = 0; // Dummy to suppress compiler warning
      printf("not recoginized");
   }

   // Pad with zeroes in beginning
   for (i = 1; i <= nd; i++)
   {
      str[i - 1] = '0'; //trim(str)
   }
   return str;
}

char *itos_unpad(int n)
{
   char str[10];
   int nd;

   // Count no. of zeroes to pad
   if (n <= 9)
   {
      std::string temp = " " + std::to_string(n);
      strcpy(str, temp.c_str());
   }
   else if (n <= 99)
   {
      std::string temp = "  " + std::to_string(n);
      strcpy(str, temp.c_str());
   }
   else if (n <= 999)
   {
      std::string temp = "   " + std::to_string(n);
      strcpy(str, temp.c_str());
   }
   else if (n <= 9999)
   {
      std::string temp = "    " + std::to_string(n);
      strcpy(str, temp.c_str());
   }
   else if (n <= 99999)
   {
      std::string temp = "     " + std::to_string(n);
      strcpy(str, temp.c_str());
   }
   else if (n <= 999999)
   {
      std::string temp = "      " + std::to_string(n);
      strcpy(str, temp.c_str());
   }
   else if (n <= 9999999)
   {
      std::string temp = "       " + std::to_string(n);
      strcpy(str, temp.c_str());
   }
   else
   {
      nd = 0; // Dummy to suppress compiler warning
      printf("not recognized unpad %d", n);
   }
}