#!/bin/bash
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

if [[ -z ${build_dir} ]]; then
  build_dir="$(mktemp -d --tmpdir build_XXXX)";
fi
if [[ -z ${src_dir} ]]; then
  if [[ -n  ${branch} ]]; then
    sudo apt-get install git;
    src_dir="$(mktemp -d --tmpdir build_XXXX)";
    git clone --branch "$1" https://github.com/statisticalbiotechnology/quandenser.git "${src_dir}/quandenser";
  else
    src_dir=$(dirname ${BASH_SOURCE})/../../../
  fi
fi
if [[ -z ${release_dir} ]]; then
  release_dir=${HOME}/release
fi

sudo apt-get update;
sudo apt-get upgrade;
sudo apt-get -y install g++ make cmake maven openjdk-8-jdk

mkdir -p ${build_dir}/tools
cd ${build_dir}/tools

if [ ! -d ${build_dir}/tools/proteowizard ]; then
  ${src_dir}/quandenser/ext/maracluster/admin/builders/install_proteowizard.sh ${build_dir}/tools
fi

mkdir -p $build_dir/quandenser
#-----cmake-----
cd $build_dir/quandenser;
echo -n "cmake quandenser.....";
cmake -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=Release -DTARBALL_BUILD=ON -DCMAKE_INSTALL_PREFIX="" -DBOOST_VARIANT_DO_NOT_USE_VARIADIC_TEMPLATES=ON -DCMAKE_PREFIX_PATH=$build_dir/tools $src_dir/quandenser;
#-----make------
echo -n "make quandenser (this will take few minutes).....";
make -j 4;
make -j 4 package;
sudo make install;

mkdir -p $release_dir
cp -v $build_dir/quandenser/quan*.tar.gz $release_dir

