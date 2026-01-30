args<-commandArgs(TRUE)
ex <- as.character(args[1])
ds<-as.character(args[2])
de<-as.character(args[3])
#print(paste(ex,ds,de))
lim_lev<-as.numeric(args[4])
max_pt <- as.numeric(args[5])
min_pt <- as.numeric(args[6])
store.folder <- as.character(args[7])
sst_anom_thresh <- as.numeric(args[8])
print(paste(ex,ds,de,lim_lev,"levels",max_pt,"-",min_pt,"km","sst anomaly",sst_anom_thresh))
# Pattern used to read data
pattern_domain <- strsplit(ex,"_")[[1]][2]
pattern_domain <- ifelse(pattern_domain=="MED","med","global")
# nsat: two data-sets available 1) allsat; 2) twosat (climate)
idx_dataset <- grep("CLIMATE",ex)
nsat <- ifelse(length(idx_dataset)>0,"twosat","allsat")

# -----------------------------
#ex<-"SWOT_WATER"
#ds<- as.character("20230801") 
#de<- as.character("20230803")
##temp_lim <- as.numeric(7) 
#lim_lev <- as.numeric(1) 
#max_pt <- as.numeric(200) 
#min_pt <- as.numeric(10) 
#store.folder <- as.character("/gpfs/work/bonaduce/NEW/ADAC/DATA/SSH_DYN/EDDY_WATER/SWOT_WATER_LOW_CY_OUT") #"./" #
#print(paste(ds,de,sep="-"))
# -----------------------------

MYD<- sprintf("/gpfs/work/bonaduce/NEW/ADAC/DATA/SSH_DYN/EDDY_WATER")
source(sprintf("%s/NAMELIST_%s.R", MYD,ex))

#NEW AC
#ssh_varname <- "ssh"
#u_varname <- "ug"
#v_varname <- "vg"
#t_varname <- "tg"
rotvel_threshold <- 0.01

.libPaths("/gpfs/home/bonaduce/R_library")
library(ncdf4,lib.loc="/gpfs/home/bonaduce/R_library")
library(SDMTools,lib.loc="/gpfs/home/bonaduce/R_library")
library(date,lib.loc="/gpfs/home/bonaduce/R_library")

starting_date<-ds # as.integer(args[2]) ## first frame read
ending_date<-de # as.integer(args[3])  ## last frame read

startdate<-as.Date(starting_date,format="%Y%m%d")
enddate<-as.Date(ending_date,format="%Y%m%d")

if(!file.exists(store.folder)) { dir.create(store.folder) }
exp.folder <- store.folder

# Check time-window
#command<-sprintf("ls -1 %s/*_%s.nc",dirONE,ex)
command<-sprintf("ls -1 %s/*_%s.nc",dirONE,ex_lab)
listfile<-system(command,inter=TRUE)
listfile<-basename(listfile)
npiece<-length(unlist(strsplit(listfile[1],split="_")))
list_days<-unlist(lapply(strsplit(listfile,split="_"), "[", 1))
list_days<-unlist(lapply(strsplit(list_days,split="\\."), "[", 1))
list_days<-as.Date(list_days,format="%Y%m%d")
idx_days <- list_days >= startdate & list_days <= enddate
list_days <- list_days[idx_days]
nd <- length(list_days)
startdate <- list_days[1]
enddate <- list_days[nd]
step_time <- as.numeric(list_days[2] - list_days[1])

#   Number of frames read
#n_times<- as.numeric(difftime(enddate,startdate,"days"))+1
n_times<- nd
#max number of eddies
#max_eddies=500
max_eddies=10000

if (ex == "ALT_GLO") {
 max_eddies <- 10000
}
pt_out <- sprintf("Max number of eddies selected : %s",max_eddies)
print(pt_out)

gc_clear<- function() { gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc()}
source(sprintf("%s/Ini.R",MYD))
source(sprintf("%s/filled_contour2.R",MYD))

        ########
        ##   output files
        #######
string_startdate<-gsub("-","",startdate)
string_enddate<-gsub("-","",enddate)

file_Ekin_edd_rel_id<-sprintf("%s/Ekin_eddy_rel_id_%s_%s.dat",exp.folder,string_startdate,string_enddate)
file_median_x<-sprintf("%s/median_x_%s_%s.dat",exp.folder,string_startdate,string_enddate)
file_median_y<-sprintf("%s/median_y_%s_%s.dat",exp.folder,string_startdate,string_enddate)
file_ncell<-sprintf("%s/numcell_%s_%s.dat",exp.folder,string_startdate,string_enddate)
file_clockwise_eddy<-sprintf("%s/clockwise_eddy_%s_%s.dat",exp.folder,string_startdate,string_enddate)
file_amplitude_edd<-sprintf("%s/amplitude_eddy_%s_%s.dat",exp.folder,string_startdate,string_enddate)
file_depth_edd<-sprintf("%s/depth_eddy_%s_%s.dat",exp.folder,string_startdate,string_enddate)
file_upwelling_edd<-sprintf("%s/upwelling_eddy_%s_%s.dat",exp.folder,string_startdate,string_enddate)

# output array
#
eddies_tot_num<-array(NA,dim=c(n_times))
Ekin_eddy_rel_id_t<-array(NA,dim=c(n_times,max_eddies))
median_x_t <-array(NA,dim=c(n_times,max_eddies))
median_y_t <-array(NA,dim=c(n_times,max_eddies))
cycl_eddy_t<-array(NA,dim=c(n_times,max_eddies))
ncell_t   <-array(NA,dim=c(n_times,max_eddies))
amplitude_eddy_t<-array(NA,dim=c(n_times,max_eddies))
depth_eddy_t<-array(NA,dim=c(n_times,max_eddies))
upwelling_eddy_t<-array(NA,dim=c(n_times,max_eddies))

        ###########################
        #   input parameter
        ###########################

# namefiles
filemask<-sprintf("%s/meshmask.nc",dirTWO)

    ## Horizontal segmentation

# max area per each eddy =100km
lim_max_pt<-max_pt*1000
lim_min_pt<-min_pt*1000
lim_lev<-lim_lev
second_set<-seq(1.5,-1.5,-0.01)
#estension for boundary treatment
ext_x_window<-1

if (ex == "ALT_GLO") {
 ext_x_window<-27
}

# VERTICAL segmentation,percentage respect to surface value
vortex_strength<- 0

#Reading dimensions, x,y extension ecc
ncm<-nc_open(filemask)
#varu_3D<-ncm$var[[u_varname]]
varu_3D<-ncm$var[["lsm"]]
varu3Dsize<-varu_3D$varsize
x_window<-varu3Dsize[1]
y_window<-varu3Dsize[2]
z_window<-1 #varu3Dsize[3]
nc_close(ncm)

        #### READING GRID, MASK
#layer by layer starting from the top
print(filemask)
ncm<-nc_open(filemask)
varu_3D<-ncm$var[["lsm"]]
#varu_3D<-ncm$var[[u_varname]]
varu3Dsize<-varu_3D$varsize
varu3Ddims<-varu_3D$ndims
print(varu3Ddims)
print("toto")
arr_st_3D<-rep(1,varu3Ddims)
arr_count_3D<-varu3Dsize

masksuv<-(ncvar_get(ncm,"lsm", start=arr_st_3D, count=arr_count_3D ))>0
#masksuv<-(ncvar_get(ncm,u_varname, start=arr_st_3D, count=arr_count_3D ))>0
depth_value<-1 #ncvar_get( ncm,"depth" )

nc_close(ncm)
masksuv<- !masksuv   #attension mask is used with 1 associated to land
#rm(tmask,umask,vmask)
gc_clear()

        ### READING DX
ncm<-nc_open(filemask)
#lon<-ncm$dim$lon$vals
#lat<-ncm$dim$latvals
depth_value<-1 #ncm$dim$depth$vals
depth_value<-depth_value[1:z_window]

nav_lon<-ncvar_get(ncm,"lon")
nav_lat<-ncvar_get(ncm,"lat")
#dx<-ncvar_get(ncm,"e1t")
#dy<-ncvar_get(ncm,"e2t")
nc_close(ncm)
dz<-1


matrix_lon<-nav_lon
matrix_lat<-nav_lat

boardcond<-matrix_lon[dim(nav_lon)[1],]
matrix_lon_up<-my_shift.up(matrix_lon,1)
matrix_lon_up[1,]<-boardcond

boardcond<-rev(matrix_lat[,dim(nav_lat)[2]]-1)
matrix_lat_left<-my_shift.left(matrix_lat,1)
matrix_lat_left[,dim(nav_lat)[2]]<-boardcond

source(sprintf("%s/hav.R",MYD))
dy<-haversine(matrix_lon,matrix_lon,matrix_lat,matrix_lat_left)*1000
dx<-haversine(matrix_lon,matrix_lon_up,matrix_lat,matrix_lat)*1000
print(dim(dy))
dx[1,]<-dx[2,]
dx[dim(dx)[1],]<-dx[dim(dx)[1]-1,]
dy[,1]<-dx[,2]
print(dim(dy)[1]-1)
#dy[,dim(dy)[2]]<-dy[,dim(dy)[1]-1]
dy[,dim(dy)[2]]<-dy[,dim(dy)[2]-1]


rm(matrix_lon,matrix_lat,matrix_lon_up,matrix_lat_left,boardcond)
gc_clear()


dx[masksuv]<-NA
dy[masksuv]<-NA
#estentinon for boundary
dx<-array(rbind(dx[(x_window-ext_x_window+1):x_window,],dx,dx[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))
dy<-array(rbind(dy[(x_window-ext_x_window+1):x_window,],dy,dy[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))

#READ MEAN DATA VALUES for anomalies - Not Used

# Array for the computation of rotational velocity

XINDEX_matr<-array(1:(ext_x_window+x_window+ext_x_window),dim=c(ext_x_window+x_window+ext_x_window,y_window))
YINDEX_matr<-t(array(1:(y_window),dim=c(y_window,ext_x_window+x_window+ext_x_window)))
x_med_matrix<-array(0,dim=c(ext_x_window+x_window+ext_x_window,y_window))
y_med_matrix<-array(0,dim=c(ext_x_window+x_window+ext_x_window,y_window))
sum_ROTVEL_matrix<-array(0,dim=c(ext_x_window+x_window+ext_x_window,y_window))
sum_ROTVEL_matrix_surf<-array(0,dim=c(ext_x_window+x_window+ext_x_window,y_window))
dataT_matrix_surf<-array(0,dim=c(ext_x_window+x_window+ext_x_window,y_window))
dataT_meanTpercentage_matrix_surf<-array(0,dim=c(ext_x_window+x_window+ext_x_window,y_window))
#
# relative frame index and eddy index
rel_time<-0
max_id_prev_edd<-0
current_date<-startdate
step_time <- 1 # days



while(current_date<=enddate)
{
string_date<-gsub("-","",current_date)
#next iteration
current_date<- as.Date(current_date + step_time,format="%Y%m%d") #next iteration

print_out<-sprintf("date %s",string_date)
print(print_out)

command<-sprintf("ls -1 %s/%s_%s.nc",dirONE,string_date,ex_lab)
fileSSH<-system(command,inter=TRUE)
fileU<-fileSSH
fileV<-fileU
#NEW AC
#fileMLD<-fileV
#ncMLD<-nc_open(fileMLD)

command<-sprintf("ls -1 %s/%s_%s.nc",dirSST,string_date,"SST_DYN")

check_exp<-length(grep("VARDYN",ex))>0

if (check_exp) {
 command<-sprintf("ls -1 %s/%s_%s.nc",dirSST,string_date,"SSH_DYN")
}


fileT<-system(command,inter=TRUE)
#fileT<-fileU
nct<-nc_open(fileT)


ncs<-nc_open(fileSSH)
vars<-ncs$var[[ssh_varname]]
varsize<-vars$varsize
vardims<-vars$ndims
arr_st<-rep(1,vardims)
arr_count<-varsize

ncu<-nc_open(fileU)
ncv<-nc_open(fileV)
nct<-nc_open(fileT)

varu_3D<-ncu$var[[u_varname]]
varu3Dsize<-varu_3D$varsize
varu3Ddims<-varu_3D$ndims
print(varu3Ddims)
print("toto")
arr_st_3D<-rep(1,varu3Ddims)
arr_count_3D<-varu3Dsize

rel_time<-rel_time+1

if (rel_time>1){start.time<-Sys.time() }

arr_st[vardims]<-1
arr_st_3D[[varu3Ddims]]<-1

#setting the mask index
arr_st_3D_dz<- arr_st_3D

datas<-ncvar_get( ncs,ssh_varname, start=arr_st, count=arr_count ) #ssh(i,j)
datau<-ncvar_get( ncu,u_varname, start=arr_st_3D, count=arr_count_3D ) #u(i+1/2,j)
datav<-ncvar_get( ncv,v_varname, start=arr_st_3D, count=arr_count_3D ) #v(i,j+1/2)
#NEW AC
datat<-ncvar_get( nct,t_varname, start=arr_st_3D, count=arr_count_3D ) #v(i,j+1/2)
#dataMLD<-ncvar_get( ncMLD,"mlotst", start=arr_st, count=arr_count ) #ssh(i,j)


#dataum<-get.var.ncdf( ncum,varid=ncum$var[["uo"]]$id, start=arr_st_3D_dz, count=arr_count_3D ) #u(i+1/2,j)
#datavm<-get.var.ncdf( ncvm,varid=ncvm$var[["vo"]]$id, start=arr_st_3D_dz, count=arr_count_3D ) #v(i,j+1/2)

masksuv<-!is.na(datau)
masksuv<- !masksuv

#datas<-datas-datasm
#datau<-datau-dataum
#datav<-datav-datavm

datas[masksuv]<-NA
datau[masksuv]<-NA
datav[masksuv]<-NA
datat[masksuv]<-NA
#dataMLD[masksuv]<-NA
#boundary
datas<-array(rbind(datas[(x_window-ext_x_window+1):x_window,],datas,datas[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))
datau<-array(rbind(datau[(x_window-ext_x_window+1):x_window,],datau,datau[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))
datav<-array(rbind(datav[(x_window-ext_x_window+1):x_window,],datav,datav[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))
#NEW AC
datat<-array(rbind(datat[(x_window-ext_x_window+1):x_window,],datat,datat[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))

datas[1,]<-NA;
datas[x_window+2*ext_x_window,]<-NA
datau[1,]<-NA;
datau[x_window+2*ext_x_window,]<-NA
datav[1,]<-NA;
datav[x_window+2*ext_x_window,]<-NA
datat[1,]<-NA;
datat[x_window+2*ext_x_window,]<-NA

#dataMLD<-array(rbind(dataMLD[(x_window-ext_x_window+1):x_window,],dataMLD,dataMLD[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))

####NEW AC
pt_out<-sprintf("datas %s %s",dim(datas)[1],dim(datas)[2])
print(pt_out)
pt_out<-sprintf("dx %s %s",dim(dx)[1],dim(dx)[2])
print(pt_out)
pt_out<-sprintf("dy %s %s",dim(dy)[1],dim(dy)[2])
print(pt_out)
pt_out<-sprintf("masksuv %s %s",dim(masksuv)[1],dim(masksuv)[2])
print(pt_out)
pt_out<-sprintf("x_window %s",x_window)
print(pt_out)
pt_out<-sprintf("ext_x_window %s" ,ext_x_window)
print(pt_out)
#
print(dim(dx[(ext_x_window+1):(x_window+ext_x_window),]))

datas<-datas-sum(datas[(ext_x_window+1):(x_window+ext_x_window),]*dx[(ext_x_window+1):(x_window+ext_x_window),]*dy[(ext_x_window+1):(x_window+ext_x_window),]*(!masksuv),na.rm=T)/sum(dx[(ext_x_window+1):(x_window+ext_x_window),]*dy[(ext_x_window+1):(x_window+ext_x_window),]*(!masksuv),na.rm=T)


#### 1st filter: velocity field checker: check the existence onf a vortex
#
#                   s1v>0 (upward is positive)
#                    △
#                    |
#                    |
#      s1u<0 <-------|------->s2u>0
#                    |
#                    |
#                    V
#                   s2v<0
#NEW AC

s2v<-my_shift.down(datav,2) #last row reversed
s1v<-my_shift.up(datav,2)   # first row reversed
s2u<-my_shift.left(datau,2)  #cyclic boundary condition
s1u<-my_shift.right(datau,2)    #cyclic boundary condition

s2v[is.na(s2v)]<-0 #otherwise conditions are not well satisfied
s2u[is.na(s2u)]<-0
s1v[is.na(s1v)]<-0
s1u[is.na(s1u)]<-0
#

Wmask1<-s2u>0 &s2v<0 & s1v>0 & s1u<0 #clockwise
Wmask2<-s2v>0 & s2u<0 & s1v<0 & s1u>0 #anticlockwise
Wmask12<- Wmask2 | Wmask1       #both


rm(s2v,s2u,s1v,s1u)
gc_clear()

# position of clockwise and anticlockwise vortexes
xy_Wmask1<-which(Wmask1==1,arr.ind=TRUE)
xy_Wmask1[,1]<-xy_Wmask1[,1]-ext_x_window
xy_Wmask2<-which(Wmask2==1,arr.ind=TRUE)
xy_Wmask2[,1]<-xy_Wmask2[,1]-ext_x_window


#
# Segmenting the horizontal extention of eddies:
# -) it starts cutting SSH from the extreme value of the module of SSH
#
#
#

# working array

ccsla_old<-array(0,c(ext_x_window+x_window+ext_x_window,y_window))
Wtot<-array(NA,c(ext_x_window+x_window+ext_x_window,y_window))
abs_datas<-abs(datas)

#
#
for( ii in second_set)
{
    Wm<-datas>ii & !is.na(datas)
    ccsla_new<-ConnCompLabel(Wm)
    idf<-factor(ccsla_new[ccsla_new>0])
    if(length(idf)==0) next
    factor_dy_id<-array(tapply(dy[ccsla_new>0],idf,mean,na.rm=TRUE))
    factor_dx_id<-array(tapply(dx[ccsla_new>0],idf,mean,na.rm=TRUE))
    max_pts_rel<-3.14*lim_max_pt/factor_dy_id*lim_max_pt/factor_dx_id
    min_pts<-3.14*lim_min_pt/factor_dy_id*lim_min_pt/factor_dx_id
    Patch<-PatchStat(ccsla_new[ccsla_new>0])
    num_cell<-Patch$n.cell
    ncelCON<-num_cell> max_pts_rel | num_cell<min_pts
    ncel_f<-Patch$patchID[ncelCON]
    Wm[(ccsla_new %in% ncel_f)]<-0
    Wm0<- Wm==0
    idf<-factor(ccsla_new[!Wm0])
    #print(length(idf))
    if(length(idf)==0) next

    id_level<-levels(idf)
    max_old_id<-array(tapply(ccsla_old[!Wm0],idf,max,na.rm=TRUE))
    min_old_id<-array(tapply(ccsla_old[!Wm0],idf,min,na.rm=TRUE))

    id_el<-id_level[!(max_old_id==min_old_id) & !is.infinite(max_old_id)]
    ccsla_new[ccsla_new %in% id_el]<-NA
    uv_f<-unique(ccsla_new[Wmask12==1 & !is.na(ccsla_new) & !Wm0 ])
    Wtot[(ccsla_new %in% uv_f)]<- 1
    ccsla_old<-ConnCompLabel(Wtot)

#image(ii)
}
revert_ssh<- -datas

for( ii in second_set)
{
    Wm<-revert_ssh>ii & !is.na(revert_ssh)

    ccsla_new<-ConnCompLabel(Wm)
    idf<-factor(ccsla_new[ccsla_new>0])
    if(length(idf)==0) next
    factor_dy_id<-array(tapply(dy[ccsla_new>0],idf,mean,na.rm=TRUE))
    factor_dx_id<-array(tapply(dx[ccsla_new>0],idf,mean,na.rm=TRUE))
    max_pts_rel<-3.14*lim_max_pt/factor_dy_id*lim_max_pt/factor_dx_id
    min_pts<-3.14*lim_min_pt/factor_dy_id*lim_min_pt/factor_dx_id
    Patch<-PatchStat(ccsla_new[ccsla_new>0])
    num_cell<-Patch$n.cell
    ncelCON<-num_cell> max_pts_rel | num_cell<min_pts
    ncel_f<-Patch$patchID[ncelCON]
    Wm[(ccsla_new %in% ncel_f)]<-0

    Wm0<- Wm==0
    idf<-factor(ccsla_new[!Wm0])
    id_level<-levels(idf)
    max_old_id<-array(tapply(ccsla_old[!Wm0],idf,max,na.rm=TRUE))
    min_old_id<-array(tapply(ccsla_old[!Wm0],idf,min,na.rm=TRUE))

    id_el<-id_level[!(max_old_id==min_old_id) & !is.infinite(max_old_id)]
    ccsla_new[ccsla_new %in% id_el]<-NA
    uv_f<-unique(ccsla_new[Wmask12==1 & !is.na(ccsla_new) & !Wm0 ])
    Wtot[(ccsla_new %in% uv_f)]<- 1
    ccsla_old<-ConnCompLabel(Wtot)

}

rm(ccsla,ccsla_old,factor_dy_id,factor_dx_id,max_pts_rel,min_pts,Patch,ncel_f,ncelCON)#
gc_clear()
########
#       final results: eddy mask is in Wtot
#       cc discriminates different eddies
########

Wtot[is.na(Wtot)]<-0
cc<-ConnCompLabel(Wtot)

Patch<-PatchStat(cc[cc>0])
shapeind<-Patch$shape.index
#constr <- median(shapeind)+(sd(shapeind))
#constr <- median(shapeind)+(2*sd(shapeind))
#constr <- median(shapeind)+(3*sd(shapeind))
#constr <- 8
constr <- 10000 # set to large value ----> consider all
ncelCON<-shapeind>constr
#ncelCON<-shapeind>8
ncel_f<-Patch$patchID[ncelCON]
Wtot[(cc %in% ncel_f)]<-0
cc<-ConnCompLabel(Wtot)

rm(Patch,ncel_f,ncelCON,constr)
gc_clear()
#image(ii*100)

# for speeding up, 0's are removed
good_cc<-!is.na(cc) & cc>0
id<-cc[good_cc]
idf<-factor(array(id))
index_ed<-as.integer(levels(idf)) # ID of each eddy


#######
#   calculating the surface properties
#####

### mean position
xy_Wtot<-which(good_cc,arr.ind=TRUE)
median_x<-array(tapply(xy_Wtot[,1],idf,mean))
median_y<-array(tapply(xy_Wtot[,2],idf,mean))

### relative velocity calculation
x_med_matrix[,]<-0
y_med_matrix[,]<-0
sum_ROTVEL_matrix_surf[,]<-0
temp<-median_x[as.integer(idf)]
x_med_matrix[as.logical(Wtot)]<- temp
temp<-median_y[as.integer(idf)]
y_med_matrix[as.logical(Wtot)]<- temp
x_matr<-round(XINDEX_matr-x_med_matrix)
y_matr<-round(YINDEX_matr-y_med_matrix)
x_matr[!as.logical(Wtot)]<-NA
y_matr[!as.logical(Wtot)]<-NA
dist_temp<-round(sqrt(x_matr^2+y_matr^2))
ROTVEL<-((datau)*y_matr*dy-x_matr*dx*(datav))/(sqrt( (y_matr*dy)^2+(x_matr*dx)^2))
sum_ROTVEL<-array(tapply(ROTVEL[good_cc],idf,mean,na.rm=TRUE))
temp<-sum_ROTVEL[as.integer(idf)]
sum_ROTVEL_matrix_surf[as.logical(Wtot)]<- temp

# Temp
dataT_matrix_surf[,]<-0
dataT_meanTpercentage_matrix_surf[,]<-0
#MLD
#max_MLD<-array(tapply(dataMLD[good_cc],idf,max,na.rm=TRUE))
#eddy depth level
#depth_eddy<-array(1,length(index_ed))

#going down to different depths, some eddies will be removed
cc_depth<-cc

cat("checking z=",1)

gc_clear()

##
## removing eddies that do not extend more than a certain value
## attention! if you use htis constrains you should re-indexing everything

#depth_eddy_in_m<-depth_value[depth_eddy]
#print(max_MLD)
#print(depth_eddy_in_m)
#depth_eddy_constrain<- depth_eddy_in_m>max_MLD
#Wtot[as.logical(Wtot)]<-depth_eddy_constrain[idf]

#Wtot[abs(sum_ROTVEL_matrix_surf)<0.01 & !is.na(sum_ROTVEL_matrix_surf)] <- 0
Wtot[abs(sum_ROTVEL_matrix_surf)<rotvel_threshold & !is.na(sum_ROTVEL_matrix_surf)] <- 0

###NEW AC
Wtot_bf_T<-Wtot

Wtot[is.na(Wtot)]<-0
cc<-ConnCompLabel(Wtot)
good_cc<-!is.na(cc) & cc>0
len_lev<-length(levels(factor(array(cc[good_cc]))))
idf_1<-factor(array(cc[good_cc]),labels=1:(len_lev))

# AB
data_to_interp<-datat
Wtemp<-Wtot
source(sprintf("%s/INTERP.R",MYD))
dataT_eddy_only<-datat-data_to_interp
sum_dataT_eddy<-array(tapply(abs(dataT_eddy_only[good_cc]),idf_1,max,na.rm=TRUE))
temp<-sum_dataT_eddy[as.integer(idf_1)]
dataT_matrix_surf[as.logical(Wtot)]<- temp

sum_meanTpercentage_eddy<-array(tapply(datat[good_cc],idf_1,max,na.rm=TRUE))
temp<-sum_meanTpercentage_eddy[as.integer(idf_1)]
dataT_meanTpercentage_matrix_surf[as.logical(Wtot)]<- temp

#
#maxT_percent<-0.05
#Wtot[abs(dataT_matrix_surf)<maxT_percent*dataT_meanTpercentage_matrix_surf & !is.na(dataT_matrix_surf)] <- 0
#Wtot[abs(dataT_matrix_surf)<0.1 & !is.na(dataT_matrix_surf)] <- 0 #  > 0.1 SST anomalies
Wtot[abs(dataT_matrix_surf) < sst_anom_thresh & !is.na(dataT_matrix_surf)] <- 0 #  > 0.1 SST anomalies

####
###

Wtot[is.na(Wtot)]<-0
cc<-ConnCompLabel(Wtot)

good_cc<-!is.na(cc) & cc>0
id<-cc[good_cc]
idf<-factor(array(id))
index_ed<-as.integer(levels(idf)) # ID of each eddy
depth_eddy<-array(1,length(index_ed))
#######
#   calculating the surface properties
#####

### mean position
xy_Wtot<-which(good_cc,arr.ind=TRUE)
median_x<-array(tapply(xy_Wtot[,1],idf,mean))
median_y<-array(tapply(xy_Wtot[,2],idf,mean))
##
if (length(depth_eddy)== 0)
{
    APP=FALSE
    #if (file.exists(file_Ekin_edd_rel_id) & rel_time==1)
    if (file.exists(file_Ekin_edd_rel_id) & rel_time>1)
    {
        APP=TRUE
    }
    write(0, file = file_Ekin_edd_rel_id,ncolumns = 1,append = APP, sep = " ")
    write(0, file = file_median_x,ncolumns = 1,append = APP, sep = " ")
    write(0, file = file_median_y,ncolumns = 1,append = APP, sep = " ")
    write(0, file = file_ncell,ncolumns = 1,append = APP, sep = " ")
    write(0, file = file_clockwise_eddy,ncolumns = 1,append = APP, sep = " ")
    write(0, file = file_amplitude_edd,ncolumns = 1,append = APP, sep = " ")
    write(0, file = file_depth_edd,ncolumns = 1,append = APP, sep = " ")
    write(0, file = file_upwelling_edd,ncolumns = 1,append = APP, sep = " ")
    rough_num_edd <- 1
    next
}


if(rel_time>1) max_id_prev_edd<-max_id_prev_edd+rough_num_edd #numbering eddies in the new frame
rough_num_edd<-length(index_ed)#not filtered

#orther properties of interests
#cycl_eddy1<-array(tapply(array(Wmask1[good_cc]),idf,sum))#clockwise
#cycl_eddy2<-array(tapply(-array(Wmask2[good_cc]),idf,sum))#counterclock
#cycl_eddy<-cycl_eddy1+cycl_eddy2
#cycl_eddy[cycl_eddy>0]<- 1 #do not switch!!
#cycl_eddy[cycl_eddy<0]<- -1

cycl_eddy<-array(tapply(array(sum_ROTVEL_matrix_surf[good_cc]),idf,sum,na.rm=TRUE))#clockwise
cycl_eddy[cycl_eddy>0]<- 1 #do not switch!!
cycl_eddy[cycl_eddy<0]<- -1

#mean SSH
mean_sla<-array(tapply(array(datas[good_cc]),idf,mean))

#eddy num cell
ncell<-PatchStat(cc[good_cc])$n.cell

#eddy amplitude
min_ssh<-array(tapply(-array(datas[good_cc]),idf,min))
max_ssh<-array(tapply(-array(datas[good_cc]),idf,max))
amplitude_eddy<-max_ssh-min_ssh

#
#point with the SSH extreme

which_extr_ssh_local<-array(tapply(array(abs_datas[good_cc]),idf,which.max))
which_extr_ssh_x<-array(0,dim=c(length(which_extr_ssh_local)))
which_extr_ssh_y<-array(0,dim=c(length(which_extr_ssh_local)))
xval_list<-array(tapply(xy_Wtot[,1],idf,array))
yval_list<-array(tapply(xy_Wtot[,2],idf,array))
for( ii in seq(1,length(which_extr_ssh_local)))
{
    temp_x<-array(unlist(xval_list[ii]))
    temp_y<-array(unlist(yval_list[ii]))
    which_extr_ssh_x[ii]<-temp_x[which_extr_ssh_local[ii]]
    which_extr_ssh_y[ii]<-temp_y[which_extr_ssh_local[ii]]
}

########
#   removing the boundary condition EW-WE
######
cond_torus<-median_x>(ext_x_window) & median_x<(ext_x_window+x_window)

index_ed<-index_ed[cond_torus]
median_x<-median_x[cond_torus]
median_y<-median_y[cond_torus]
ncell<-ncell[cond_torus]
Ekin_eddy_rel_id<-index_ed
cycl_eddy<-cycl_eddy[cond_torus]
amplitude_eddy<-amplitude_eddy[cond_torus]
depth_eddy<-depth_eddy[cond_torus]
mean_sla<-mean_sla[cond_torus]
which_extr_ssh_x<-which_extr_ssh_x[cond_torus]
which_extr_ssh_y<-which_extr_ssh_y[cond_torus]
which_extr_ssh<-array(cbind(which_extr_ssh_x,which_extr_ssh_y),dim=c(length(which_extr_ssh_x),2))

sla_on_extreme<-datas[which_extr_ssh]
upwelling_eddy<-array(-1,length(sla_on_extreme))
upwelling_eddy[sla_on_extreme>mean_sla]<- 1

rm(Wmask12)
gc_clear()

#print image
#image(rel_time)

########### end of the single frame ########
##
#           STORING VALUES
###
if(length(Ekin_eddy_rel_id)>max_eddies)
   {
      print("number of eddies larger than max...increasing")
        max_eddies_old<- max_eddies
        max_eddies<-length(Ekin_eddy_rel_id)

        Ekin_eddy_rel_id_t_old<-Ekin_eddy_rel_id_t
        Ekin_eddy_rel_id_t<-array(NA,dim=c(n_times,max_eddies))
        Ekin_eddy_rel_id_t[1:rel_time,1:max_eddies_old]<-Ekin_eddy_rel_id_t_old[1:rel_time,1:max_eddies_old]
        rm(Ekin_eddy_rel_id_t_old)

        median_x_t_old<-median_x_t
        median_x_t<-array(NA,dim=c(n_times,max_eddies))
        median_x_t[1:rel_time,1:max_eddies_old]<-median_x_t_old[1:rel_time,1:max_eddies_old]
        rm(median_x_t_old)

        median_y_t_old<-median_y_t
        median_y_t<-array(NA,dim=c(n_times,max_eddies))
        median_y_t[1:rel_time,1:max_eddies_old]<-median_y_t_old[1:rel_time,1:max_eddies_old]
        rm(median_y_t_old)

        cycl_eddy_t_old<-cycl_eddy_t
        cycl_eddy_t<-array(NA,dim=c(n_times,max_eddies))
        cycl_eddy_t[1:rel_time,1:max_eddies_old]<-cycl_eddy_t_old[1:rel_time,1:max_eddies_old]
        rm(cycl_eddy_t_old)

        ncell_t_old<-ncell_t
        ncell_t<-array(NA,dim=c(n_times,max_eddies))
        ncell_t[1:rel_time,1:max_eddies_old]<-ncell_t_old[1:rel_time,1:max_eddies_old]
        rm(ncell_t_old)

        amplitude_eddy_t_old<-amplitude_eddy_t
        amplitude_eddy_t<-array(NA,dim=c(n_times,max_eddies))
        amplitude_eddy_t[1:rel_time,1:max_eddies_old]<-amplitude_eddy_t_old[1:rel_time,1:max_eddies_old]
        rm(amplitude_eddy_t_old)

	depth_eddy_t_old<-depth_eddy_t
        depth_eddy_t<-array(NA,dim=c(n_times,max_eddies))
        depth_eddy_t[1:rel_time,1:max_eddies_old]<-depth_t_old[1:rel_time,1:max_eddies_old]
        rm(depth_eddy_t_old)

        upwelling_eddy_t_old<-upwelling_eddy_t
        upwelling_eddy_t<-array(NA,dim=c(n_times,max_eddies))
        upwelling_eddy_t[1:rel_time,1:max_eddies_old]<-upwelling_eddy_t_old[1:rel_time,1:max_eddies_old]
        rm(upwelling_eddy_t_old)


    }
#
#output arrays with all the values
#
eddies_tot_num[rel_time]<-length(Ekin_eddy_rel_id)
Ekin_eddy_rel_id_t[rel_time,1:eddies_tot_num[rel_time]]<-Ekin_eddy_rel_id
median_x_t[rel_time,1:eddies_tot_num[rel_time]]<-median_x-ext_x_window #positions are now the right one
median_y_t[rel_time,1:eddies_tot_num[rel_time]]<-median_y
cycl_eddy_t[rel_time,1:eddies_tot_num[rel_time]]<-cycl_eddy
ncell_t[rel_time,1:eddies_tot_num[rel_time]]<-ncell
amplitude_eddy_t[rel_time,1:eddies_tot_num[rel_time]]<-amplitude_eddy
depth_eddy_t[rel_time,1:eddies_tot_num[rel_time]]<-depth_eddy
upwelling_eddy_t[rel_time,1:eddies_tot_num[rel_time]]<-upwelling_eddy

# creating the files
    if (rel_time==1)
        {
        write(eddies_tot_num[rel_time], file = file_Ekin_edd_rel_id,ncolumns = 1, sep = " ")
        write(Ekin_eddy_rel_id_t[rel_time,1:eddies_tot_num[rel_time]], file = file_Ekin_edd_rel_id,ncolumns = length(Ekin_eddy_rel_id),append = TRUE, sep = " ")
        write(eddies_tot_num[rel_time], file = file_median_x,ncolumns = 1, sep = " ")
        write(median_x_t[rel_time,1:eddies_tot_num[rel_time]], file = file_median_x,ncolumns = length(median_x),append = TRUE, sep = " ")
        write(eddies_tot_num[rel_time], file = file_median_y,ncolumns = 1, sep = " ")
        write(median_y_t[rel_time,1:eddies_tot_num[rel_time]], file = file_median_y,ncolumns = length(median_y),append = TRUE, sep = " ")
        write(eddies_tot_num[rel_time], file = file_ncell,ncolumns = 1, sep = " ")
        write(ncell_t[rel_time,1:eddies_tot_num[rel_time]], file = file_ncell,ncolumns = length(ncell),append = TRUE, sep = " ")
        write(eddies_tot_num[rel_time], file = file_clockwise_eddy,ncolumns = 1, sep = " ")
        write(cycl_eddy_t[rel_time,1:eddies_tot_num[rel_time]], file = file_clockwise_eddy,ncolumns = length(cycl_eddy),append = TRUE, sep = " ")
        write(eddies_tot_num[rel_time], file = file_amplitude_edd,ncolumns = 1, sep = " ")
        write(amplitude_eddy_t[rel_time,1:eddies_tot_num[rel_time]], file = file_amplitude_edd,ncolumns = length(amplitude_eddy),append = TRUE, sep = " ")
        write(eddies_tot_num[rel_time], file = file_depth_edd,ncolumns = 1, sep = " ")
        write(depth_eddy_t[rel_time,1:eddies_tot_num[rel_time]], file = file_depth_edd,ncolumns = length(depth_eddy),append = TRUE, sep = " ")
        write(eddies_tot_num[rel_time], file = file_upwelling_edd,ncolumns = 1, sep = " ")
        write(upwelling_eddy_t[rel_time,1:eddies_tot_num[rel_time]], file = file_upwelling_edd,ncolumns = length(upwelling_eddy),append = TRUE, sep = " ")

# --------
# Write NC
# --------
dimX <- ncdim_def("x", "", 1:x_window,create_dimvar=FALSE)
dimY <- ncdim_def("y", "", 1:y_window,create_dimvar=FALSE)

Wtot_excluded_ind<- ncvar_def("Wtot_exc", "counts", list(dimX,dimY), 9.96921e+36, longname=" eddy removed by T analysis ", prec="single",compression=4)
Wtot_final_ind<- ncvar_def("Wtot", "counts", list(dimX,dimY), 9.96921e+36, longname=" eddy retained ", prec="single",compression=4)
ssh_ind<- ncvar_def("SSH", "m", list(dimX,dimY), 9.96921e+36, longname=" ssh ", prec="single",compression=4)
temp_ind<- ncvar_def("temp", "degC", list(dimX,dimY), 9.96921e+36, longname=" temperture ", prec="single",compression=4)
temp_noeddy_ind<- ncvar_def("temp_noeddy", "degC", list(dimX,dimY), 9.96921e+36, longname=" temp_noeddy ", prec="single",compression=4)
T_anomaly_excluded_ind<- ncvar_def("T_exc_anom", "degC", list(dimX,dimY), 9.96921e+36, longname=" excluded eddy temperature anomaly ", prec="single",compression=4)
T_anomaly_ind<- ncvar_def("T_anom", "degC", list(dimX,dimY), 9.96921e+36, longname=" eddy temperature anomaly ", prec="single",compression=4)
T_anomaly_cy_ind<- ncvar_def("T_anom_cy", "degC", list(dimX,dimY), 9.96921e+36, longname=" eddy temperature anomaly ", prec="single",compression=4)
T_anomaly_anticy_ind<- ncvar_def("T_anom_anticy", "degC", list(dimX,dimY), 9.96921e+36, longname=" eddy temperature anomaly ", prec="single",compression=4)
clock_ind<-ncvar_def("clockwise", "adim", list(dimX,dimY), 9.96921e+36, longname=" clockwise/counter ", prec="single",compression=4)
nav_lon_ind<-ncvar_def("nav_lon", "degrees_east", list(dimX,dimY), 9.96921e+36, longname="", prec="single",compression=4)
nav_lat_ind<-ncvar_def("nav_lat", "degrees_north", list(dimX,dimY), 9.96921e+36, longname=" ", prec="single",compression=4)

file_out<-sprintf("%s/Detected_eddies_%s.nc",exp.folder,string_date)
ncout<-nc_create( file_out, list(ssh_ind,temp_ind,temp_noeddy_ind,Wtot_excluded_ind,Wtot_final_ind,T_anomaly_excluded_ind,T_anomaly_ind,T_anomaly_cy_ind,T_anomaly_anticy_ind,clock_ind,nav_lat_ind,nav_lon_ind),force_v4=TRUE )
start_arr<-c(1,1)
count_arr<-c(x_window,y_window)
ncvar_put( ncout, ssh_ind, datas[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, temp_ind, datat[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, temp_noeddy_ind, data_to_interp[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, nav_lat_ind, nav_lat, start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, nav_lon_ind, nav_lon, start=start_arr, count=count_arr, verbose=FALSE )
ttemp<-Wtot_bf_T-Wtot
ncvar_put( ncout, Wtot_excluded_ind, ttemp[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, Wtot_final_ind, Wtot[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )

ttemp<-dataT_eddy_only*Wtot
#ttemp<-T_anomaly_excluded_ind*Wtot
ncvar_put( ncout, T_anomaly_ind, ttemp[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ttemp<-dataT_eddy_only*Wtot*(sum_ROTVEL_matrix_surf>0)
ncvar_put( ncout, T_anomaly_cy_ind, ttemp[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ttemp<-dataT_eddy_only*Wtot*(sum_ROTVEL_matrix_surf<0)
ncvar_put( ncout, T_anomaly_anticy_ind, ttemp[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
#ttemp<-T_anomaly_excluded_ind*(1-Wtot)
ttemp<-dataT_eddy_only*(1-Wtot)
ncvar_put( ncout, T_anomaly_excluded_ind, ttemp[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ttemp<-sum_ROTVEL_matrix_surf*Wtot
ttemp[ttemp>0]<- 1
ttemp[ttemp<0]<- -1
ncvar_put( ncout, clock_ind, ttemp[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )

nc_close(ncout)


# ----------
        next
        }
        cat("rel time=",rel_time)

#open files are write frame by frame

write(eddies_tot_num[rel_time], file = file_Ekin_edd_rel_id,ncolumns = 1,append = TRUE, sep = " ")
write(Ekin_eddy_rel_id_t[rel_time,1:eddies_tot_num[rel_time]], file = file_Ekin_edd_rel_id,ncolumns = length(Ekin_eddy_rel_id),append = TRUE, sep = " ")

write(eddies_tot_num[rel_time], file = file_median_x,ncolumns = 1,append = TRUE, sep = " ")
write(median_x_t[rel_time,1:eddies_tot_num[rel_time]], file = file_median_x,ncolumns = length(median_x),append = TRUE, sep = " ")

write(eddies_tot_num[rel_time], file = file_median_y,ncolumns = 1,append = TRUE, sep = " ")
write(median_y_t[rel_time,1:eddies_tot_num[rel_time]], file = file_median_y,ncolumns = length(median_y),append = TRUE, sep = " ")

write(eddies_tot_num[rel_time], file = file_ncell,ncolumns = 1,append = TRUE, sep = " ")
write(ncell_t[rel_time,1:eddies_tot_num[rel_time]], file = file_ncell,ncolumns = length(ncell),append = TRUE, sep = " ")

write(eddies_tot_num[rel_time], file = file_clockwise_eddy,ncolumns = 1,append = TRUE, sep = " ")
write(cycl_eddy_t[rel_time,1:eddies_tot_num[rel_time]], file = file_clockwise_eddy,ncolumns = length(cycl_eddy),append = TRUE, sep = " ")

write(eddies_tot_num[rel_time], file = file_amplitude_edd,ncolumns = 1,append = TRUE, sep = " ")
write(amplitude_eddy_t[rel_time,1:eddies_tot_num[rel_time]], file = file_amplitude_edd,ncolumns = length(amplitude_eddy),append = TRUE, sep = " ")

write(eddies_tot_num[rel_time], file = file_depth_edd,ncolumns = 1,append = TRUE, sep = " ")
write(depth_eddy_t[rel_time,1:eddies_tot_num[rel_time]], file = file_depth_edd,ncolumns = length(depth_eddy),append = TRUE, sep = " ")

write(eddies_tot_num[rel_time], file = file_upwelling_edd,ncolumns = 1,append = TRUE, sep = " ")
write(upwelling_eddy_t[rel_time,1:eddies_tot_num[rel_time]], file = file_upwelling_edd,ncolumns = length(upwelling_eddy),append = TRUE, sep = " ")

gc_clear()

end.time<-Sys.time()
exec.time<-end.time-start.time
print(exec.time)

nc_close(ncu)
nc_close(ncv)
nc_close(ncs)
nc_close(nct)

#NEW AC
dimX <- ncdim_def("x", "", 1:x_window,create_dimvar=FALSE)
dimY <- ncdim_def("y", "", 1:y_window,create_dimvar=FALSE)

Wtot_excluded_ind<- ncvar_def("Wtot_exc", "counts", list(dimX,dimY), 9.96921e+36, longname=" eddy removed by T analysis ", prec="single",compression=4)
Wtot_final_ind<- ncvar_def("Wtot", "counts", list(dimX,dimY), 9.96921e+36, longname=" eddy retained ", prec="single",compression=4)
ssh_ind<- ncvar_def("SSH", "m", list(dimX,dimY), 9.96921e+36, longname=" ssh ", prec="single",compression=4)
temp_ind<- ncvar_def("temp", "degC", list(dimX,dimY), 9.96921e+36, longname=" temperture ", prec="single",compression=4)
temp_noeddy_ind<- ncvar_def("temp_noeddy", "degC", list(dimX,dimY), 9.96921e+36, longname=" temp_noeddy ", prec="single",compression=4)
T_anomaly_excluded_ind<- ncvar_def("T_exc_anom", "degC", list(dimX,dimY), 9.96921e+36, longname=" excluded eddy temperature anomaly ", prec="single",compression=4)
T_anomaly_ind<- ncvar_def("T_anom", "degC", list(dimX,dimY), 9.96921e+36, longname=" eddy temperature anomaly ", prec="single",compression=4)
T_anomaly_cy_ind<- ncvar_def("T_anom_cy", "degC", list(dimX,dimY), 9.96921e+36, longname=" eddy temperature anomaly ", prec="single",compression=4)
T_anomaly_anticy_ind<- ncvar_def("T_anom_anticy", "degC", list(dimX,dimY), 9.96921e+36, longname=" eddy temperature anomaly ", prec="single",compression=4)
clock_ind<-ncvar_def("clockwise", "adim", list(dimX,dimY), 9.96921e+36, longname=" clockwise/counter ", prec="single",compression=4)
nav_lon_ind<-ncvar_def("nav_lon", "degrees_east", list(dimX,dimY), 9.96921e+36, longname="", prec="single",compression=4)
nav_lat_ind<-ncvar_def("nav_lat", "degrees_north", list(dimX,dimY), 9.96921e+36, longname=" ", prec="single",compression=4)

file_out<-sprintf("%s/Detected_eddies_%s.nc",exp.folder,string_date)
ncout<-nc_create( file_out, list(ssh_ind,temp_ind,temp_noeddy_ind,Wtot_excluded_ind,Wtot_final_ind,T_anomaly_excluded_ind,T_anomaly_ind,T_anomaly_cy_ind,T_anomaly_anticy_ind,clock_ind,nav_lat_ind,nav_lon_ind),force_v4=TRUE )
start_arr<-c(1,1)
count_arr<-c(x_window,y_window)
ncvar_put( ncout, ssh_ind, datas[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, temp_ind, datat[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, temp_noeddy_ind, data_to_interp[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, nav_lat_ind, nav_lat, start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, nav_lon_ind, nav_lon, start=start_arr, count=count_arr, verbose=FALSE )
ttemp<-Wtot_bf_T-Wtot
ncvar_put( ncout, Wtot_excluded_ind, ttemp[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, Wtot_final_ind, Wtot[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )

ttemp<-dataT_eddy_only*Wtot
#ttemp<-T_anomaly_excluded_ind*Wtot
ncvar_put( ncout, T_anomaly_ind, ttemp[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ttemp<-dataT_eddy_only*Wtot*(sum_ROTVEL_matrix_surf>0)
ncvar_put( ncout, T_anomaly_cy_ind, ttemp[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ttemp<-dataT_eddy_only*Wtot*(sum_ROTVEL_matrix_surf<0)
ncvar_put( ncout, T_anomaly_anticy_ind, ttemp[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
#ttemp<-T_anomaly_excluded_ind*(1-Wtot)
ttemp<-dataT_eddy_only*(1-Wtot)
ncvar_put( ncout, T_anomaly_excluded_ind, ttemp[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ttemp<-sum_ROTVEL_matrix_surf*Wtot
ttemp[ttemp>0]<- 1
ttemp[ttemp<0]<- -1
ncvar_put( ncout, clock_ind, ttemp[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )

nc_close(ncout)

}


nc_close(ncm)
#close.ncdf(ncum)
#close.ncdf(ncvm)
