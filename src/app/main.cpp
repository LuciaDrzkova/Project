// Isosurface Singularity Lab - entry point of the interactive viewer.

#include "IsosurfaceViewer.h"

#include <exception>
#include <iostream>

int main()
{
    try
    {
        iso_app::IsosurfaceViewer viewer("Isosurface Singularity Lab", 1280, 800);
        return viewer.run();
    }
    catch (const std::exception& error)
    {
        std::cerr << "Isosurface Singularity Lab failed: " << error.what() << '\n';
        return 1;
    }
}
