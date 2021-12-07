::: Centralized place for urls and files for all windows builders ...
::: please do not change compression type in urls, since decompression is
::: hardcoded in the respective buiding scripts


::: maven
set MVN_BASE=apache-maven-3.6.3
set MVN_URL=http://ftp.unicamp.br/pub/apache/maven/maven-3/3.6.3/binaries/%MVN_BASE%-bin.zip

::: java
set JAVA_BASE=openjdk-8u212-b04
set JAVA_URL=https://github.com/AdoptOpenJDK/openjdk8-upstream-binaries/releases/download/jdk8u212-b04/OpenJDK8u-x64_windows_8u212b04.zip


EXIT /B
