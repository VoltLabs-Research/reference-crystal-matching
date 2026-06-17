from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps, cmake_layout


class ReferenceCrystalMatchingConan(ConanFile):
    name = "reference-crystal-matching"
    version = "1.0.1"
    package_type = "static-library"
    license = "MIT"
    settings = "os", "arch", "compiler", "build_type"
    default_options = {"hwloc/*:shared": True}
    requires = (
        "boost/1.88.0",
        "onetbb/2021.12.0",
        "coretoolkit/[>=2.5]",
        "structure-identification/[>=2.1]",
        "spdlog/1.14.1",
        "nlohmann_json/3.11.3",
        "yaml-cpp/0.8.0",
    )
    exports_sources = "CMakeLists.txt", "include/*", "src/*", "dependencies/*", "structures/*"

    def layout(self):
        cmake_layout(self)

    def generate(self):
        CMakeToolchain(self).generate()
        CMakeDeps(self).generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.set_property(
            "cmake_target_name",
            "reference-crystal-matching::reference-crystal-matching",
        )
        self.cpp_info.libs = ["reference-crystal-matching_lib"]
        self.cpp_info.requires = [
            "boost::headers",
            "onetbb::onetbb",
            "coretoolkit::coretoolkit",
            "structure-identification::structure-identification",
            "nlohmann_json::nlohmann_json",
            "spdlog::spdlog",
            "yaml-cpp::yaml-cpp",
        ]
