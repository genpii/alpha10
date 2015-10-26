#pragma once
#include <windows.h>

HINSTANCE hInstance = GetModuleHandle(0);
HWND hwnd;
HDC hdc;
PAINTSTRUCT ps;
HPEN hpen;
HBRUSH brush;
int x_coordinate=0, y_coordinate=0;

int get_x(), get_y();

LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp) {

	switch (msg) {
	default:
		break;
	case WM_DESTROY:
		PostQuitMessage(0);
		return 0;
	case WM_MOUSEMOVE:
		x_coordinate = get_x();
		y_coordinate = get_y();
		break;
	}
	return DefWindowProc(hwnd, msg, wp, lp);
}

LRESULT CALLBACK btn(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp) {

	switch (msg) {
	case WM_LBUTTONUP:
		return 1;
	default:
		return 0;
	}
}

int pageb(int x1, int y1, int nx, int ny, char *ch) {

	WNDCLASS winc;

	winc.style = CS_HREDRAW | CS_VREDRAW;
	winc.lpfnWndProc = WndProc;
	winc.cbClsExtra = winc.cbWndExtra = 0;
	winc.hInstance = hInstance;
	winc.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	winc.hCursor = LoadCursor(NULL, IDC_ARROW);
	winc.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);
	winc.lpszMenuName = NULL;
	winc.lpszClassName = TEXT("window");

	if (!RegisterClass(&winc)) return 0;

	hwnd = CreateWindow(TEXT("window"), ch, WS_OVERLAPPEDWINDOW|WS_VISIBLE, x1, y1, nx, ny, NULL, NULL, hInstance, NULL);

	if (hwnd == NULL) return 0;
	
	return 1;
}

void gsline(int x1, int y1, int x2, int y2, COLORREF color) {

	hpen = CreatePen(PS_SOLID, 0, color);
	SelectObject(hdc, hpen);
	MoveToEx(hdc, x1, y1, NULL);
	LineTo(hdc, x2, y2);
	DeleteObject(hpen);
}

void gsline_d(int x1, int y1, int x2, int y2, COLORREF color) {

	hpen = CreatePen(PS_DASH, 0, color);
	SelectObject(hdc, hpen);
	MoveToEx(hdc, x1, y1, NULL);
	LineTo(hdc, x2, y2);
	DeleteObject(hpen);
}

void gsrect(int x1, int y1, int x2, int y2, COLORREF color_in, COLORREF color_out) {

	hpen = CreatePen(PS_SOLID, 0, color_out);
	SelectObject(hdc, hpen);
	brush = CreateSolidBrush(color_in);
	SelectObject(hdc, brush);
	Rectangle(hdc, x1, y1, x2, y2);
	DeleteObject(hpen);
	DeleteObject(brush);
}

COLORREF gscol256(int ir, int ig, int ib) {

	COLORREF color;

	color = RGB(ir, ig, ib);
	return color;
}

void ptext(int x1, int y1, int len, char *ch) {

	TextOut(hdc, x1, y1, ch, len);
}

int get_x() {

	POINT pt;
	GetCursorPos(&pt);
	ScreenToClient(hwnd, &pt);
	return pt.x;
}

int get_y() {

	POINT pt;
	GetCursorPos(&pt);
	ScreenToClient(hwnd, &pt);
	return pt.y;
}

void cursorpos(int x, int y) {

	POINT pt;
	pt.x = x;
	pt.y = y;
	ClientToScreen(hwnd, &pt);
	SetCursorPos(pt.x, pt.y);
}

void eraseg() {

	EndPaint(hwnd, &ps);
	InvalidateRect(hwnd, NULL, TRUE);
	hdc = BeginPaint(hwnd, &ps);
}

int messageloop(void) {

	WNDCLASS winc;
	MSG msg;

	//winc.lpfnWndProc = WndProc;
	while (GetMessage(&msg, NULL, 0, 0)>0) {
		DispatchMessage(&msg);
	}
	return msg.wParam;
}

void start(void) {

	hdc = BeginPaint(hwnd, &ps);
}

void end(void) {

	EndPaint(hwnd, &ps);
}

//gscolor
//0:black, 1:white, 2:blue, 3:red