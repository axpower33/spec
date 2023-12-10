// Spec.cpp : Определяет точку входа для приложения.
//

#include "framework.h"
#include "array.h"
#include "Spec.h"
#include <stdlib.h>
#include <string.h>
#include <conio.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <dos.h>
#include <iostream>
using namespace std;

typedef unsigned int uint;

float** As;
float** P;
float** U;
float* W;
float* Xcl;
float* Ycl;
float* Zcl;
float* X;
float* L;
float* del;
float* XX;
float* YY;
float* ZZ;
float* Chi;
float* ya;
float* db;
int countCall;

const float M_PI = 3.14;
const int N = 50;
const int nn = 3 * N;
const int Kmx = 49;
const float Xmn = -5e17;
const float Xmx = 5e17;
float eh = 1.78;
float Rm = 12e-7;
float c = 3e10;
float hb = 1.05e-34;
float wp = 2 * M_PI * c / 1.361e-5;
float gb = 0.0019 * wp;
float vf = 0.0047 * c;
float gRm = gb + vf / Rm;
float Lmn = 2e-5;
float Lmx = 1e-4;
float Lm = 4.1e-5;

struct complex
{
    float Re, Im;
} e[Kmx], et[Kmx];

float w[Kmx], lev[Kmx], n[Kmx], k[Kmx];

float sqr(float x) noexcept
{
    return x * x;
}

bool FormArray() 
{
    bool mem = false;

    if (!(As = (float**)makeAr(nn, (nn + 1) / 2 + 1, sizeof(float)))) mem = true;


    if (!(P = (float**)makeAr(nn, nn, sizeof(float)))) mem = true;


    if (!(W = new float[nn])) mem = true;


    if (!(U = (float**)makeAr(nn, nn, sizeof(float)))) mem = true;


    if (!(Xcl = new float[N])) mem = true;

    if (!(Ycl = new float[N])) mem = true;


    if (!(Zcl = new float[N])) mem = true;


    if (!(XX = new float[Kmx])) mem = true;


    if (!(YY = new float[Kmx])) mem = true;


    if (!(ZZ = new float[Kmx])) mem = true;


    if (!(Chi = new float[Kmx])) mem = true;


    if (!(X = new float[Kmx])) mem = true;


    if (!(L = new float[Kmx])) mem = true;


    if (!(del = new float[Kmx])) mem = true;


    if (!(ya = new float[Kmx])) mem = true;


    if (!(db = new float[Kmx])) mem = true;

    if (!mem) return true;

    return false;
}

void DestrAr() noexcept
{
    delete[] Xcl;
    delete[] Ycl;
    delete[] Zcl;
    delete[] W;
    delAr((void**)U);
    delAr((void**)As);
    delAr((void**)P);
    delete[] XX;
    delete[] YY;
    delete[] ZZ;
    delete[] Chi;
    delete[] X;
}

float funMax(float* fun)
{
    float Cm;
    Cm = fun[0];
    for (int k = 1; k < Kmx; k++) if (Cm < fun[k]) Cm = fun[k];
    return Cm;
}

float funMin(float* fun)
{
    float Cm;
    Cm = fun[0];
    for (int k = 1; k < Kmx; k++) if (Cm > fun[k]) Cm = fun[k];
    return Cm;
}

void GraphIm(HDC hdc, float* X, float* fun) 
{
    float RNG = 0.17;
    _flushall();


    for (int k = 1; k < Kmx; k++) {
        MoveToEx(hdc, X[k - 1], fun[k - 1], NULL);
        LineTo(hdc, X[k], fun[k]);
    }
    _getch();

}

void Deigen(int n, int mv) noexcept
{
    const float range = 0.00001;
    uint iq, mm, ll, lq, mq, jq, ind, ij, i, j;
    uint l, m, lm, ilq, imq, im, il, imr, ilr, lmm;
    float anorm, x, y, sinx, sinx2, cosx, cosx2, sincs;
    float thr, anrmx;
    if (mv != 1)
    {
        iq = -n;
        for (j = 1; j < n + 1; j++)
        {
            iq += n;
            for (uint i = 1; i < n + 1; i++)
                if (i == j) P[(iq + i - 1) / nn][(iq + i - 1) % nn] = 1;
                else P[(iq + i - 1) / nn][(iq + i - 1) % nn] = 0;
        }
    }
    anorm = 0.0;
    for (j = 2; j < n + 1; j++)
        for (i = 1; i < j; i++)
        {
            ij = i + (j * j - j) / 2;
            anorm += pow(As[(ij - 1) / nn][(ij - 1) % nn], 2);
        }
    if (anorm > 0)
    {
        anorm = sqrt(anorm) * 1.414;
        anrmx = anorm * range / n;
        // compute threshold:THR
        ind = 0; thr = anorm;
        // Їp®ўҐpЄ  ­®p¬л ¬ вpЁжл
        do
        {
            thr /= n;

        m50: l = 1;
        m55: m = l + 1;
            //   compute SIN and COS

        m60: mq = (m * m - m) / 2;
            lq = (l * l - l) / 2;
            lm = l + mq;
            if (fabs(As[(lm - 1) / nn][(lm - 1) % nn]) >= thr)
            {
                ind = 1;
                ll = l + lq;
                mm = m + mq;
                x = 0.5 * (As[(ll - 1) / nn][(ll - 1) % nn] - As[(mm - 1) / nn][(mm - 1) % nn]);
                y = -As[(lm - 1) / nn][(lm - 1) % nn] / sqrt(pow(As[(lm - 1) / nn][(lm - 1) % nn], 2) + x * x);
                if (x < 0) y = -y;
                sinx = y / sqrt(2.0 * (1.0 + (sqrt(1.0 - y * y))));
                sinx2 = sinx * sinx;
                cosx = sqrt(1.0 - sinx2);
                cosx2 = cosx * cosx;
                sincs = sinx * cosx;
                //     rotate L and M columns
                ilq = n * (l - 1);
                imq = n * (m - 1);
                for (i = 1; i < n + 1; i++)
                {
                    iq = (i * i - i) / 2;
                    if (i != l)
                    {
                        if (i < m) { im = i + mq; goto m95; }
                        else { if (i == m) goto m115; im = m + iq; }
                    m95:     if (i < l) { il = i + lq; goto m110; }
                        il = l + iq;
                    m110:    x = As[(il - 1) / nn][(il - 1) % nn] * cosx - As[(im - 1) / nn][(im - 1) % nn] * sinx;
                        As[(im - 1) / nn][(im - 1) % nn] = As[(il - 1) / nn][(il - 1) % nn] * sinx + As[(im - 1) / nn][(im - 1) % nn] * cosx;
                        As[(il - 1) / nn][(il - 1) % nn] = x;
                    }
                m115:  if (mv != 1)
                {
                    ilr = ilq + i;
                    imr = imq + i;
                    x = P[(ilr - 1) / nn][(ilr - 1) % nn] * cosx - P[(imr - 1) / nn][(imr - 1) % nn] * sinx;
                    P[(imr - 1) / nn][(imr - 1) % nn] = P[(ilr - 1) / nn][(ilr - 1) % nn] * sinx + P[(imr - 1) / nn][(imr - 1) % nn] * cosx;
                    P[(ilr - 1) / nn][(ilr - 1) % nn] = x;
                }
                }
                x = 2.0 * As[(lm - 1) / nn][(lm - 1) % nn] * sincs;
                y = As[(ll - 1) / nn][(ll - 1) % nn] * cosx2 + As[(mm - 1) / nn][(mm - 1) % nn] * sinx2 - x;
                x = As[(ll - 1) / nn][(ll - 1) % nn] * sinx2 + As[(mm - 1) / nn][(mm - 1) % nn] * cosx2 + x;
                As[(lm - 1) / nn][(lm - 1) % nn] = (As[(ll - 1) / nn][(ll - 1) % nn] - As[(mm - 1) / nn][(mm - 1) % nn]) * sincs + As[(lm - 1) / nn][(lm - 1) % nn] * (cosx2 - sinx2);
                As[(ll - 1) / nn][(ll - 1) % nn] = y;
                As[(mm - 1) / nn][(mm - 1) % nn] = x;
            }
            //   test for completion
            if (m != n) { m++; goto m60; }
            if (l != (n - 1)) { l++; goto m55; }
            if (ind == 1) { ind = 0; goto m50; }
        } while (thr > anrmx);
    }
    for (i = 1; i < nn + 1; i++)
    {
        int k = (i * (i + 1)) / 2;
        W[i - 1] = As[(k - 1) / nn][(k - 1) % nn];
    }
    m = 0; l = 1;
    for (i = 1; i < nn * nn + 1; i++)
    {
        m++;
        if (m == nn + 1)
        {
            m = 1; l++;
        }
        U[m - 1][l - 1] = P[(i - 1) / nn][(i - 1) % nn];
    }


}

void tranc()
{
    
    FILE* cl;
    fopen_s(&cl, "sp.dat", "r");

    for (int i = 0; i < Kmx; i++)
    {
        fscanf_s(cl, "%e %e %e ", lev[i], n[i], k[i]);

        w[i] = lev[i] * 1.6e-19 / hb;
        L[i] = 2 * M_PI * c / w[i];

        et[i].Re = sqr(n[i]) - sqr(k[i]);
        et[i].Im = 2 * n[i] * k[i];

        e[i].Re = et[i].Re + sqr(wp) *
            (1 / (sqr(w[i]) + sqr(gb)) - 1 / (sqr(w[i]) + sqr(gRm)));
        e[i].Im = et[i].Im + sqr(wp) / w[i] *
            (-gb / (sqr(w[i]) + sqr(gb)) + gRm / (sqr(w[i]) + sqr(gRm)));

        ya[i] = 1.0 / (Rm * Rm * Rm) * 3 * eh * e[i].Im / (sqr(e[i].Re - eh) + sqr(e[i].Im));
        db[i] = 2.0 / 3 * pow(2 * M_PI / L[i], 3);
        del[i] = ya[i] + db[i];

        X[i] = -1.0 / (Rm * Rm * Rm) * (1 + 3 * eh * (e[i].Re - eh) / (sqr(e[i].Re - eh) + sqr(e[i].Im)));
    }
    fclose(cl);
}

void mainris(HDC hdc) 
{
    uint i, j;
    float z, h;
    float Ux, Uy, Uz;
    float Uux, Uuy, Uuz;


    char fdat[128];
    fdat[0] = 0;
    HRGN hrgn = CreateRectRgn(0, 0, 640, 480);
    HBRUSH hbr = CreateSolidBrush(RGB(0, 0, 0));
    FillRgn(hdc, hrgn, hbr);
    
    FILE* F1;
    fopen_s((FILE**)&F1,"frsp.dat", "rt");
    for (i = 0; i < N; i++)
            fscanf_s(F1, "%e %e %e", &Xcl[i], &Ycl[i], &Zcl[i]);
        
    fclose(F1);
    
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            if (i != j)
            {
                float rx = Xcl[i] - Xcl[j];
                float ry = Ycl[i] - Ycl[j];
                float rz = Zcl[i] - Zcl[j];
                float r2 = rx * rx + ry * ry + rz * rz;
                float r12 = sqrt(r2) / 1.65;
                float r5 = 1.0 / (r2 * r2 * r12);
                U[3 * i][3 * j] = (r2 - 3 * rx * rx) * r5;
                U[3 * i][3 * j + 1] = -3 * rx * ry * r5;
                U[3 * i][3 * j + 2] = -3 * rx * rz * r5;
                U[3 * i + 1][3 * j] = U[3 * i][3 * j + 1];
                U[3 * i + 1][3 * j + 1] = (r2 - 3 * ry * ry) * r5;
                U[3 * i + 1][3 * j + 2] = -3 * ry * rz * r5;
                U[3 * i + 2][3 * j] = U[3 * i][3 * j + 2];
                U[3 * i + 2][3 * j + 1] = U[3 * i + 1][3 * j + 2];
                U[3 * i + 2][3 * j + 2] = (r2 - 3 * rz * rz) * r5;
            }
            else
            {
                U[3 * i][3 * j] = 0;
                U[3 * i][3 * j + 1] = 0;
                U[3 * i][3 * j + 2] = 0;
                U[3 * i + 1][3 * j] = 0;
                U[3 * i + 1][3 * j + 1] = 0;
                U[3 * i + 1][3 * j + 2] = 0;
                U[3 * i + 2][3 * j] = 0;
                U[3 * i + 2][3 * j + 1] = 0;
                U[3 * i + 2][3 * j + 2] = 0;
            }
    int kl = 0;
    for (i = 0; i < nn; i++)
        for (j = 0; j < nn; j++)
            if (i >= j)
            {
                kl++;
                As[(kl - 1) / nn][(kl - 1) % nn] = U[i][j];
                //   printf("\n<%u|W|%u>=%e",(i+1)/2,(j+1)/2,U[i][j]);
            }
    Deigen(nn, 0);
    int k;
    for (k = 0; k < Kmx; k++)
    {
        XX[k] = 0;
        YY[k] = 0;
        ZZ[k] = 0;
    }
    //   clrscr();
    tranc();
    for (k = 0; k < Kmx; k++)
    { //X[k]=Xmn+k*(Xmx-Xmn)/Kmx;
  
        for (int n = 0; n < nn; n++)
        {
            z = -X[k] + W[n];
            for (i = 0; i < N; i++)
            {
                Uux = U[3 * i][n];
                Uuy = U[3 * i + 1][n];
                Uuz = U[3 * i + 2][n];
                h = z * z + del[k] * del[k];
                for (j = 0; j < N; j++)
                {
                    Ux = Uux * U[3 * j][n] / h;
                    Uy = Uuy * U[3 * j + 1][n] / h;
                    Uz = Uuz * U[3 * j + 2][n] / h;
                    XX[k] = XX[k] + del[k] * Ux;
                    YY[k] = YY[k] + del[k] * Uy;
                    ZZ[k] = ZZ[k] + del[k] * Uz;
                }
            }
        }
    }

    for (k = 0; k < Kmx; k++)
    {
        XX[k] = XX[k] / (3 * N);
        YY[k] = YY[k] / (3 * N);
        ZZ[k] = ZZ[k] / (3 * N);
        Chi[k] = (XX[k] + YY[k] + ZZ[k]) * 4 * 2 * M_PI / (L[k] * Rm * Rm);
    }
    GraphIm(hdc, L, Chi);
    DestrAr();
}
#define MAX_LOADSTRING 100

// Глобальные переменные:
HINSTANCE hInst;                                // текущий экземпляр
WCHAR szTitle[MAX_LOADSTRING];                  // Текст строки заголовка
WCHAR szWindowClass[MAX_LOADSTRING];            // имя класса главного окна

// Отправить объявления функций, включенных в этот модуль кода:
ATOM                MyRegisterClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);

int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPWSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    // TODO: Разместите код здесь.

    // Инициализация глобальных строк
    LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
    LoadStringW(hInstance, IDC_SPEC, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);

    // Выполнить инициализацию приложения:
    if (!InitInstance (hInstance, nCmdShow))
    {
        return FALSE;
    }

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_SPEC));

    MSG msg;

    // Цикл основного сообщения:
    while (GetMessage(&msg, nullptr, 0, 0))
    {
        if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }

    return (int) msg.wParam;
}



//
//  ФУНКЦИЯ: MyRegisterClass()
//
//  ЦЕЛЬ: Регистрирует класс окна.
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style          = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc    = WndProc;
    wcex.cbClsExtra     = 0;
    wcex.cbWndExtra     = 0;
    wcex.hInstance      = hInstance;
    wcex.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_SPEC));
    wcex.hCursor        = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground  = (HBRUSH)(COLOR_WINDOW+1);
    wcex.lpszMenuName   = MAKEINTRESOURCEW(IDC_SPEC);
    wcex.lpszClassName  = szWindowClass;
    wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}

//
//   ФУНКЦИЯ: InitInstance(HINSTANCE, int)
//
//   ЦЕЛЬ: Сохраняет маркер экземпляра и создает главное окно
//
//   КОММЕНТАРИИ:
//
//        В этой функции маркер экземпляра сохраняется в глобальной переменной, а также
//        создается и выводится главное окно программы.
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
   hInst = hInstance; // Сохранить маркер экземпляра в глобальной переменной

   HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
      CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);

   if (!hWnd)
   {
      return FALSE;
   }

   ShowWindow(hWnd, nCmdShow);
   UpdateWindow(hWnd);

   return TRUE;
}

//
//  ФУНКЦИЯ: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  ЦЕЛЬ: Обрабатывает сообщения в главном окне.
//
//  WM_COMMAND  - обработать меню приложения
//  WM_PAINT    - Отрисовка главного окна
//  WM_DESTROY  - отправить сообщение о выходе и вернуться
//
//
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    switch (message)
    {
    case WM_COMMAND:
        {
            int wmId = LOWORD(wParam);
            // Разобрать выбор в меню:
            switch (wmId)
            {
            case IDM_ABOUT:
                DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
                break;
            case IDM_EXIT:
                DestroyWindow(hWnd);
                break;
            default:
                return DefWindowProc(hWnd, message, wParam, lParam);
            }
        }
        break;
    case WM_PAINT:
        {
            PAINTSTRUCT ps;
            HDC hdc = BeginPaint(hWnd, &ps);
            mainris(hdc);
            // TODO: Добавьте сюда любой код прорисовки, использующий HDC...
            EndPaint(hWnd, &ps);
        }
        break;
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }
    return 0;
}

// Обработчик сообщений для окна "О программе".
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
        {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}
