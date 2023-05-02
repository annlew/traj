# TrajectoryScripts
This repository includes trajectory scripts for computing and plotting forecast trajectories 

The atos branch is prepared for running on Atos hpc/ecs and requires access to the lagranto installation of Moritz Pickl and Oden location from Linus Manusson's ftp server. Upload of data is done to box.

## Preparations

### Lagranto installation

Add the following to ~/.profile or ~/.bashrc
```
# Lagranto settings
export model="ecmwf"
export mode="ive"
export DYN_TOOLS=/home/nemp/programs/eth_tools
export LAGRANTO=${DYN_TOOLS}/lagranto.${model}-${mode}
export PATH=$PATH:${DYN_TOOLS}/bin
export PATH=$PATH:${LAGRANTO}/bin
export PATH=$PATH:/home/nemp/programs/bin

LD_LIBRARY_PATH=/usr/local/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/home/nemp/programs/eth_tools/local/lib"
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/usr/local/apps/netcdf4/4.7.4/GNU/8.4/lib"

export LD_LIBRARY_PATH
```

### Fetch Oden location

Requires access to Linus ftp server. Add the following in ~/.netrc

```
machine bolftp.ecmwf.int login lmagnusson password XXXXXX
```

### Set up box

IN TERMINAL ON ATOS
```
$ rclone config
Current remotes:

Name                 Type
====                 ====

e) Edit existing remote
n) New remote
d) Delete remote
r) Rename remote
c) Copy remote
s) Set configuration password
q) Quit config
e/n/d/r/c/s/q> n
name> box
Type of storage to configure.
Enter a string value. Press Enter for the default ("").
Choose a number from below, or type in your own value
 1 / 1Fichier
   \ "fichier"
 2 / Alias for an existing remote
   \ "alias"
 3 / Amazon Drive
   \ "amazon cloud drive"
 4 / Amazon S3 Compliant Storage Provider (AWS, Alibaba, Ceph, Digital Ocean, Dreamhost, IBM COS, Minio, Tencent COS, etc)
   \ "s3"
 5 / Backblaze B2
   \ "b2"
 6 / Box
   \ "box"
 7 / Cache a remote
   \ "cache"
 8 / Citrix Sharefile
   \ "sharefile"
 ........
Storage> box
** See help for box backend at: https://rclone.org/box/ **

OAuth Client Id
Leave blank normally.
Enter a string value. Press Enter for the default ("").
client_id>
OAuth Client Secret
Leave blank normally.
Enter a string value. Press Enter for the default ("").
client_secret>
Box App config.json location
Leave blank normally.

Leading `~` will be expanded in the file name as will environment variables such as `${RCLONE_CONFIG_DIR}`.

Enter a string value. Press Enter for the default ("").
box_config_file>
Box App Primary Access Token
Leave blank normally.
Enter a string value. Press Enter for the default ("").
access_token>

Enter a string value. Press Enter for the default ("user").
Choose a number from below, or type in your own value
 1 / Rclone should act on behalf of a user
   \ "user"
 2 / Rclone should act on behalf of a service account
   \ "enterprise"
box_sub_type>
Edit advanced config? (y/n)
y) Yes
n) No (default)
y/n> n
Remote config
Use auto config?
 * Say Y if not sure
 * Say N if you are working on a remote or headless machine
y) Yes (default)
n) No
y/n> n
For this to work, you will need rclone available on a machine that has
a web browser available.

For more help and alternate methods see: https://rclone.org/remote_setup/

Execute the following on the machine with the web browser (same rclone
version recommended):

    rclone authorize "box"

Then paste the result below:
result>
```

NEW TERMINAL ON LOCAL MACHINE

```
$ rclone authorize "box"
If your browser doesn't open automatically go to the following link: http://127.0.0.1:53682/auth?state=gnpuF2MifOxCzmI1lXogug
Log in and authorize rclone for access
Waiting for code...
Got code
Paste the following into your remote machine --->
{"access_token":"xxxx","token_type":"bearer","refresh_token":"xxxx","expiry":"2023-04-28T14:09:18.143003794+02:00"}
<---End paste
```

IN TERMINAL ON ATOS
```
Then paste the result below:
result> {"access_token":"xXXXXXXXXXXXXXXX","token_type":"bearer","refresh_token":"XXXXXXXXXXXXX","expiry":"2023-04-28T14:09:18.143003794+02:00"}
--------------------
[annas_box]
token = {"access_token":"XXXXXXXXXXXXXXXXXX","token_type":"bearer","refresh_token":"XXXXXXXXXXXXXXXXXXxx","expiry":"2023-04-28T14:09:18.143003794+02:00"}
--------------------
y) Yes this is OK (default)
e) Edit this remote
d) Delete this remote
y/e/d> y
Current remotes:

Name                 Type
====                 ====
box                  box

e) Edit existing remote
n) New remote
d) Delete remote
r) Rename remote
c) Copy remote
s) Set configuration password
q) Quit config
e/n/d/r/c/s/q> q
```


## Running in interactive mode







