# Usage

```
cd <pandora_folder>/portable_binary_builder
./build_portable_binary.sh
```

The portable binary will be in `<pandora_folder>/pandora-linux-precompiled`


# Troubleshooting

If the `phusion/holy-build-box-64:2.0.1` docker version does not work, e.g. with error:

```
+ yum install wget git -y
Loaded plugins: fastestmirror, ovl
Setting up Install Process
Determining fastest mirrors
YumRepo Error: All mirror URLs are not using ftp, http[s] or file.
 Eg. Invalid release/repo/arch combination/
removing mirrorlist with no valid mirrors: /var/cache/yum/x86_64/6/base/mirrorlist.txt
Error: Cannot find a valid baseurl for repo: base
```

build the container yourself:

```
git clone https://github.com/phusion/holy-build-box && cd holy-build-box && sudo make 64
```
