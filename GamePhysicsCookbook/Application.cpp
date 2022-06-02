#include "Application.h"

#include <raylib.h>

#include "Config.h"

Application::Application()
{
	windowWidth = config.GetIntValue(PROGRAM_CATEGORY, "width");
	windowHeight = config.GetIntValue(PROGRAM_CATEGORY, "height");
}

void Application::Run()
{
	InitWindow(windowWidth, windowHeight, config.GetTextValue(PROGRAM_CATEGORY, "name"));

	Start();

	while (!WindowShouldClose())
	{
		Update(GetFrameTime());

		BeginDrawing();
		ClearBackground(RAYWHITE);
		Draw();
		EndDrawing();
	}

	OnDestroy();

	CloseWindow();
}

void Application::Start()
{
	// Runs before the first update loop of the application
}

void Application::Update(float _dt)
{

}

void Application::Draw()
{

}

void Application::OnDestroy()
{
	// Runs after the update loop has been exited
}
