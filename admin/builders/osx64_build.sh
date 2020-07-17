#!/bin/bash
# Requirements are:
# XCode
# Command line tools (check if installed with "xcode-select -p", otherwise install with "xcode-select --install")
# MacPorts or homebrew as package manager
# PackageMaker (https://developer.apple.com/downloads ::search::
#----------------------------------------

# managing input arguments
while getopts “s:b:r:t:” OPTION; do
  case $OPTION in
    s) src_dir=${OPTARG};;
    t) branch=${OPTARG};;
    r) release_dir=${OPTARG};;
    b) build_dir=${OPTARG};;
    \?) echo "Invalid option: -${OPTARG}" >&2;;
  esac
done

if [[ ! -d /Applications/XCode.app ]]
  then
    echo "Apple developer tools are required (Search for XCode in the App Store)"
    exit 1
fi

package_manager_installed=true
if [[ -d /opt/local/var/macports ]]
  then
    echo "[ Package manager ] : MacPorts "
    package_manager="sudo port"
    other_packages="cmake gnutar wget maven32 openjdk8"
elif [[ -f ${HOME}/bin/brew ]] || [[ -f /usr/local/bin/brew ]]
  then
    echo "[ Package manager ] : Homebrew "
    if [[ -f ${HOME}/bin/brew ]]
      then
        package_manager=$HOME/bin/brew
    else
      package_manager=brew
    fi
    other_packages="cmake gnu-tar wget libomp maven"
    ${package_manager} update || true # brew.rb raises an error on the vagrant box, just ignore it
    ${package_manager} tap AdoptOpenJDK/openjdk
    ${package_manager} cask install adoptopenjdk8
else
    package_manager_installed=false
fi

if [ "$package_manager_installed" == false ]
  then
  echo "Error: no suitable package manager installed"
  echo " Get homebrew or macports:"
  echo "  Homebrew: http://brew.sh/ "
  echo "  MacPorts: http://www.macports.org/install.php"
  exit 1
fi

if [[ -z ${build_dir} ]]; then
  build_dir="$(mktemp -d -t build)";
fi
if [[ -z ${src_dir} ]]; then
  if [[ -n  ${branch} ]]; then
    if [[ ! -f /usr/bin/git ]]; then
      $package_manager install git;
    fi
    src_dir="$(mktemp -d -t src)";
    git clone --branch "$1" https://github.com/statisticalbiotechnology/quandenser.git "${src_dir}/quandenser";
  else
    # Might not work if we have symlinks in the way
    src_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../" && pwd )
  fi
fi
if [[ -z ${release_dir} ]]; then
  release_dir=${HOME}/release
fi


echo "The Builder $0 is building the Quandenser packages with src=${src_dir} an\
d build=${build_dir} for user" `whoami`
$package_manager install $other_packages

#----------------------------------------

mkdir -p ${build_dir}/tools
cd ${build_dir}/tools

if [ ! -d ${build_dir}/tools/proteowizard ]; then
  ${src_dir}/quandenser/ext/maracluster/admin/builders/install_proteowizard.sh ${build_dir}/tools
fi

# brew has trouble installing maven+openjdk for old XCode versions
if mvn -v; then
  echo "Maven installed successfully"
else
  echo "Could not find maven, installing standalone version"
  wget --quiet https://mirror.softaculous.com/apache/maven/maven-3/3.6.3/binaries/apache-maven-3.6.3-bin.zip
  unzip apache-maven-3.6.3-bin.zip
  export PATH="$PATH:${build_dir}/tools/apache-maven-3.6.3/bin"
fi

#-------------------------------------------

mkdir -p $build_dir/quandenser
#-----cmake-----
# we need to install to /usr/local instead of /usr: https://github.com/Benjamin-Dobell/Heimdall/issues/291
cd $build_dir/quandenser;
echo -n "cmake quandenser.....";
cmake -DCMAKE_CXX_COMPILER="/usr/bin/clang++" -DTARGET_ARCH="x86_64" -DCMAKE_BUILD_TYPE=Release -DBoost_COMPILER=-xgcc42 -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_PREFIX_PATH="${build_dir}/tools/" $src_dir/quandenser;
#-----make------
echo -n "make quandenser (this will take few minutes).....";
make -j 4;
make -j 4 package;
#sudo make install;

#--------------------------------------------

echo "build directory was : ${build_dir}";

mkdir -p $release_dir
cp -v $build_dir/quandenser/quan*.pkg $release_dir

