-- Nori Ray Tracer xmake configuration
-- Converted from CMakeLists.txt with reduced ext/ dependencies

set_project("nori")
set_version("1.0.0")
set_languages("c++17")

-- Add build modes
add_rules("mode.debug", "mode.release")

-- Set default to release
if not has_config("mode") then
    set_config("mode", "release")
end

-- Enable compile commands export for IDEs
-- set_policy("build.export_compile_commands", true)

-- Add main include directory
add_includedirs("include")

-- System package dependencies
add_requires("eigen", "pugixml", "openexr", "glfw", "nanogui","stb")

-- For NanoGUI and NanoVG, we might need to handle them separately
-- Option 1: Use system packages if available
-- add_requires("nanogui", "nanovg")  -- if available in package manager

-- Option 2: Build from source in ext/ (if must keep)
-- We'll assume header-only libraries remain in ext/
add_includedirs("ext/tinyformat")      -- Header-only formatting library
add_includedirs("ext/pcg32")           -- Header-only random number generator
add_includedirs("ext/hypothesis")      -- Header-only hypothesis testing
add_includedirs("ext/filesystem")      -- Header-only filesystem API

add_includedirs("ext/tbb")             -- Header-only image write (if STB is header-only)

-- If NanoGUI/NanoVG must be built from ext/, add them here
-- add_includedirs("ext/nanogui/include")
-- add_includedirs("ext/nanovg/src")

-- Platform-specific settings
if is_plat("linux") then
    add_syslinks("GL", "GLU", "X11", "Xxf86vm", "Xrandr", "pthread", "Xi", "Xinerama", "Xcursor", "dl")
elseif is_plat("macosx") then
    add_frameworks("OpenGL", "Cocoa", "IOKit", "CoreVideo")
elseif is_plat("windows") then
    add_syslinks("opengl32", "glu32", "user32", "gdi32", "shell32")
    add_defines("NOMINMAX", "_CRT_SECURE_NO_WARNINGS")
end

-- Compiler settings
if is_plat("linux", "macosx") then
    add_cxflags("-Wall", "-Wextra", "-Wno-unused-parameter")
    if is_mode("debug") then
        add_cxflags("-g", "-O0")
        add_defines("DEBUG")
    else
        add_cxflags("-O3", "-march=native")
        add_defines("NDEBUG")
    end
    
    -- Clang-specific flags
    if is_plat("linux") then
        add_cxflags("-Wno-gnu-anonymous-struct", "-Wno-c99-extensions", 
                   "-Wno-nested-anon-types", "-Wno-deprecated-register")
    end
        add_cxflags("-fcolor-diagnostics")
        add_cxflags("-fdiagnostics-color=always")
end

-- If NanoGUI needs to be built from ext/
-- target("nanogui")
--     set_kind("static")
--     add_files("ext/nanogui/src/*.cpp")
--     add_headerfiles("ext/nanogui/include/nanogui/*.h")
--     add_packages("glfw", "glew", "eigen3")

-- Main executable: nori
target("nori")
    set_kind("binary")
    
    -- Add system packages
    add_packages("eigen", "pugixml", "openexr", "glfw",  "glew","stb","nanogui")

    -- If NanoGUI is built as a target
    -- add_deps("nanogui")
    
    -- Header files (for IDE support)
    add_headerfiles("include/nori/*.h")
    
    -- Core source files
    add_files("src/bitmap.cpp")
    add_files("src/block.cpp") 
    add_files("src/accel.cpp")
    add_files("src/chi2test.cpp")
    add_files("src/common.cpp")
    add_files("src/diffuse.cpp")
    add_files("src/gui.cpp")
    add_files("src/independent.cpp")
    add_files("src/main.cpp")
    add_files("src/mesh.cpp")
    add_files("src/obj.cpp")
    add_files("src/object.cpp")
    add_files("src/parser.cpp")
    add_files("src/perspective.cpp")
    add_files("src/proplist.cpp")
    add_files("src/rfilter.cpp")
    add_files("src/scene.cpp")
    add_files("src/ttest.cpp")
    add_files("src/warp.cpp")
    add_files("src/microfacet.cpp")
    add_files("src/mirror.cpp")
    add_files("src/dielectric.cpp")
    
    -- Custom added integrators
    add_files("src/normal.cpp")
    add_files("src/simple.cpp") 
    add_files("src/ao.cpp")
    add_files("src/area.cpp")
    add_files("src/whitted.cpp")
    add_files("src/path_mats.cpp")
    add_files("src/path_ems.cpp")
    
    -- Windows-specific: link with zlib
    if is_plat("windows") then
        add_packages("zlib")
    end

-- Test executable: warptest
target("warptest")
    set_kind("binary")

    add_packages("eigen", "pugixml", "openexr", "glfw",  "glew","stb","nanogui")

    -- If NanoGUI is built as a target
    -- add_deps("nanogui")
    
    add_files("src/warp.cpp")
    add_files("src/warptest.cpp")
    add_files("src/microfacet.cpp")
    add_files("src/object.cpp")
    add_files("src/proplist.cpp")
    add_files("src/common.cpp")
    add_files("src/bitmap.cpp")

-- Custom tasks
task("setup")
    set_menu({
        description = "Setup minimal required header-only libraries from ext/",
        usage = "xmake setup"
    })
    on_run(function ()
        print("Setting up minimal header-only libraries...")
        local header_only_libs = {
            "tinyformat",  -- String formatting
            "pcg32",       -- Random number generator
            "hypothesis",  -- Statistical tests
            "filesystem",  -- Filesystem API
            "stb"          -- STB image libraries
        }
        
        for _, lib in ipairs(header_only_libs) do
            if not os.isdir("ext/" .. lib) then
                print("Warning: Missing ext/" .. lib)
                print("Please ensure header-only library '" .. lib .. "' is in ext/")
            end
        end
        
        print("\nFor GUI support, you may need to:")
        print("1. Install NanoGUI from system packages, or")
        print("2. Build NanoGUI from ext/nanogui (uncomment relevant sections in xmake.lua)")
        print("\nNow run: xmake config && xmake")
    end)

task("install-deps")
    set_menu({
        description = "Install system dependencies",
        usage = "xmake install-deps"
    })
    on_run(function ()
        if is_plat("linux") then
            print("Installing dependencies on Linux...")
            print("Please run ONE of the following:")
            print("\nUbuntu/Debian:")
            print("  sudo apt install libeigen3-dev libpugixml-dev libopenexr-dev \\")
            print("                   libglfw3-dev libglew-dev libgl1-mesa-dev \\")
            print("                   libglu1-mesa-dev libx11-dev libtbb-dev")
            print("\nFedora:")
            print("  sudo dnf install eigen3-devel pugixml-devel openexr-devel \\")
            print("                   glfw-devel glew-devel mesa-libGL-devel \\")
            print("                   mesa-libGLU-devel libX11-devel tbb-devel")
            print("\nArch:")
            print("  sudo pacman -S eigen pugixml openexr glfw-x11 glew \\")
            print("                 mesa glu libx11 intel-tbb")
            print("\nFor NanoGUI (optional):")
            print("  Check if available in your distribution or build from source")
        elseif is_plat("macosx") then
            print("Installing dependencies on macOS...")
            print("Please run:")
            print("  brew install eigen pugixml openexr glfw glew tbb")
            print("\nFor NanoGUI (optional):")
            print("  brew install nanogui  # if available")
        elseif is_plat("windows") then
            print("On Windows, xmake will try to install packages automatically.")
            print("You may need to install:")
            print("  - Visual Studio with C++ support")
            print("  - vcpkg or other package managers for dependencies")
        end
    end)

-- Optional: Build configuration summary
task("info")
    set_menu({
        description = "Show build configuration",
        usage = "xmake info"
    })
    on_run(function ()
        print("Nori Ray Tracer Build Configuration")
        print("===================================")
        print("Platform: " .. os.host())
        print("Architecture: " .. os.arch())
        print("Build mode: " .. (is_mode("debug") and "debug" or "release"))
        print("\nRequired packages:")
        print("  - eigen3: Linear algebra")
        print("  - pugixml: XML parsing") 
        print("  - openexr: HDR image support")
        print("  - glfw: Window management")
        print("  - glew: OpenGL extensions")
        print("  - tbb: Threading support")
        print("\nHeader-only libraries (in ext/):")
        print("  - tinyformat: String formatting")
        print("  - pcg32: Random numbers")
        print("  - hypothesis: Statistical tests")
        print("  - filesystem: File operations")
        print("  - stb: Image I/O")
        print("\nOptional:")
        print("  - nanogui: User interface")
        print("  - nanovg: Vector graphics")
    end)
