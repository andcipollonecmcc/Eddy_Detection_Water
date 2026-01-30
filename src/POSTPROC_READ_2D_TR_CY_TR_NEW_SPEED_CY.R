args<-commandArgs(TRUE)
ex <- as.character(args[1]) #"BLKSEA" #
ds<-as.character(args[2]) # "20080101" #
de<-as.character(args[3]) # "20171231" #
ds_proc<- as.character(args[4]) #"20101010" #
de_proc<- as.character(args[5])# "20101014" #
temp_lim <- as.numeric(args[6]) #3 #
lim_lev <- as.numeric(args[7]) #23 #
max_pt <- as.numeric(args[8]) #200 #
min_pt <- as.numeric(args[9]) #5 #
store.folder <- as.character(args[10]) #"./" #
print(paste(ds,de,"lifetime",temp_lim,"days"))
# Pattern used to read data
pattern_domain <- strsplit(ex,"_")[[1]][2]
pattern_domain <- ifelse(pattern_domain=="MED","med","global")
# nsat: two data-sets available 1) allsat; 2) twosat (climate)
idx_dataset <- grep("CLIMATE",ex)
nsat <- ifelse(length(idx_dataset)>0,"twosat","allsat")

# -----------------------------
#ex<-"SWOT_WATER"
#ds<-as.character("20230727")
#de<-as.character("20240501") 
#ds_proc<- as.character("20230727") 
#de_proc<- as.character("20230803")
#temp_lim <- as.numeric(7) 
#lim_lev <- as.numeric(1) 
#max_pt <- as.numeric(200) 
#min_pt <- as.numeric(10) 
#store.folder <- as.character("/gpfs/work/bonaduce/NEW/ADAC/DATA/SSH_DYN/EDDY_WATER/SWOT_WATER_OUT") #"./" #
#print(paste(ds,de,"lifetime",temp_lim,"days"))

#ds<-as.character("20230801")
#de<-as.character("20230803") 
#ds_proc<- as.character("20230802") 
#de_proc<- as.character("20230803")
#temp_lim <- as.numeric(2) 
#lim_lev <- as.numeric(1) 
#max_pt <- as.numeric(200) 
#min_pt <- as.numeric(10) 
#store.folder <- as.character("/gpfs/work/bonaduce/NEW/ADAC/DATA/SSH_DYN/EDDY_WATER/SWOT_WATER_LOW_CY_OUT")
#print(paste(ds,de,"lifetime",temp_lim,"days"))


# -----------------------------

MYD<-sprintf("/gpfs/work/bonaduce/NEW/ADAC/DATA/SSH_DYN/EDDY_WATER")
namelist<-sprintf("%s/NAMELIST_%s.R",MYD,ex)
print(namelist)
source(namelist)

#NEW AC
#ssh_varname <- "ssh"
#u_varname <- "ug"
#v_varname <- "vg"
#t_varname <- "tg"
rotvel_threshold <- 0.01


.libPaths("/gpfs/home/bonaduce/R_library")
library(ncdf4,lib.loc="/gpfs/home/bonaduce/R_library")
library(SDMTools,lib.loc="/gpfs/home/bonaduce/R_library")
#library(fields,lib.loc="/gpfs/home/bonaduce/R_library")
library(date,lib.loc="/gpfs/home/bonaduce/R_library")
suppressMessages(library(matrixcalc,lib.loc="/gpfs/home/bonaduce/R_library"))



string_startdate<-ds_proc
string_enddate<-de_proc

starting_date_experiment<-ds
ending_date_experiment<-de
starting_date<-ds_proc # as.integer(args[2]) ## first frame read
ending_date<-de_proc # as.integer(args[3])  ## last frame read

startdate_experiment<-as.Date(starting_date_experiment,format="%Y%m%d")
enddate_experiment<-as.Date(ending_date_experiment,format="%Y%m%d")
startdate<-as.Date(starting_date,format="%Y%m%d")
enddate<-as.Date(ending_date,format="%Y%m%d")

seqDATE<-seq(startdate,enddate, by="day")

dirMEAN=dirTWO

MYD<- sprintf("/gpfs/work/bonaduce/NEW/ADAC/DATA/SSH_DYN/EDDY_WATER")
if(!file.exists(store.folder)) { dir.create(store.folder) }
exp.folder <- store.folder
#exp.folder <- sprintf("%s_POS",exp.folder)
print(exp.folder)
if(!file.exists(exp.folder)) { dir.create(exp.folder) }
###########################
#   input parameter
###########################

## Horizontal segmentation

# max area per each eddy =100km
lim_max_pt<-max_pt*1000
lim_min_pt<-min_pt*1000

# ssh is cut from 0.5 to 0.1 according to following resolution
#second_set<-c(c(0.7,0.6,0.5,0.45,0.40,0.30),seq(0.25,-0.25,-0.02),c(-0.3,-0.4,-0.45,-0.5,-0.6,-0.7))
#second_set<-seq(3.3,-3.3,-0.01)
second_set<-seq(1.5,-1.5,-0.01)
#estension for boundary treatment
ext_x_window<-1
if (ex == "ALT_GLO") {
 ext_x_window<-27
}
# VERTICAL segmentation,percentage respect to surface value
#vortex_strength<- 0.25
vortex_strength<- 0
#

# Check time-window
# Check time-window in Arctic/Antarctic Data
#command<-sprintf("ls -1 %s/dt_*_multimission_sea_level_*.nc",dirONE)
#command<-sprintf("ls -1 /gpfs/work/bonaduce/AVISO_ARCTIC_CREG025_UPD/SPLIT/dt_*_multimission_sea_level_*.nc")

#startdate <- as.Date("20181230",format="%Y%m%d") # list_days[1]
#enddate <- as.Date("20181231",format="%Y%m%d") #list_days[nd]
startdate <- as.Date(ds_proc,format="%Y%m%d") # list_days[1]
enddate <- as.Date(de_proc,format="%Y%m%d") #list_days[nd]
step_time <- 1 #as.numeric(list_days[2] - list_days[1])  # dt= 3 days

#
# converting date to number of frame
#
	#
	#  current frame
	#framread_first<-as.numeric(difftime(startdate,startdate_experiment,"days"))+1
	#framread_last<-as.numeric(difftime(enddate,startdate_experiment,"days"))+1
	n_times_rel<- as.numeric(difftime(enddate,startdate,"days"))+1
	# total frame considered
	#n_times<-as.numeric(difftime(enddate_experiment,startdate_experiment,"days"))+1

	# Update Variables: framread_first, framread_last,n_times,ntime_rel
	#framread_first <- which(startdate==list_days_all)
	#framread_last <- which(enddate==list_days_all)
	#list_days_experiment<-seq(as.Date("20140101",format="%Y%m%d"),as.Date("20181231",format="%Y%m%d"),by="day")
	list_days_experiment<-seq(as.Date(ds,format="%Y%m%d"),as.Date(de,format="%Y%m%d"),by="day")
    framread_first <- which(startdate==list_days_experiment)
	framread_last <- which(enddate==list_days_experiment)
	# Update n_time
	#n_times_rel <- nd
	n_times <- length(list_days_experiment)
	#max number of eddies
	#max_eddies=500
	max_eddies=5000

	if (ex == "ALT_GLO") {
	 max_eddies <- 10000
	 max_eddies <- 1000
	 max_eddies <- 500
	}
	pt_out <- sprintf("Max number of eddies selected : %s",max_eddies)
	print(pt_out)

	#####
	#temp_lim<-2 #days
	days_each_record<-1
	#days_each_record<-3 # Arctic data: available every 3 days

	###
	###  FOR EDDY PROFILE
	lifetime<-temp_lim
	###

	gc_clear<- function() { gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc()}
	source(sprintf("%s/Ini.R",MYD))
	source(sprintf("%s/filled_contour2.R",MYD))
	##

	###
	###
	###

	# namefiles
	filemask<-sprintf("%s/meshmask.nc",dirTWO)

	#Reading dimensions, x,y extension ecc
	ncm<-nc_open(filemask)
	#varu_3D<-ncm$var[[u_varname]]
	varu_3D<-ncm$var[["lsm"]]
	varu3Dsize<-varu_3D$varsize
	x_window<-varu3Dsize[1]
	y_window<-varu3Dsize[2]
	#z_window<-1 #varu3Dsize[3]
	nc_close(ncm)

	z_window <- lim_lev

#print("Read Masks")
# Output File
#fo<-sprintf("%s/MASK_REGIONS.nc",dirTWO)

#nc<-nc_open(fo)

#nv<-vector()
#for (n in 1:nc$nvars) {
# nv[n]<-nc$var[[n]]$name
#}
#idx_var<-grep("MASK",nv)
#nv<-nv[idx_var]
#num_reg<-length(nv)
#
#MASK_REGION<-array(NA,dim=c(x_window+2*ext_x_window,y_window,length(nv)))
#
#for (nn in 1:num_reg) {
# MASK_REGION[(ext_x_window+1):(x_window+ext_x_window),,nn ]<-ncvar_get(nc,nv[nn])
## #E-W bounday treatments
# MASK_REGION[1:ext_x_window,,nn ]<-MASK_REGION[(x_window+1):(x_window+ext_x_window),,nn]
# MASK_REGION[(x_window+ext_x_window+1):(x_window+2*ext_x_window),,nn ]<-MASK_REGION[(ext_x_window+1):(2*ext_x_window),,nn ]
#}
#MASK_REGION[is.na(MASK_REGION)]<-0
#MASK_REGION[MASK_REGION != 0 ]<-1
#nc_close(nc)

		########
		##   output files
		#######

	file_Ekin_edd_rel_id<-sprintf("%s/Ekin_eddy_rel_id.dat",exp.folder)
	file_eddy_index<-sprintf("%s/eddy_index.dat",exp.folder)
	file_median_x<-sprintf("%s/median_x.dat",exp.folder)
	file_median_y<-sprintf("%s/median_y.dat",exp.folder)
	file_ncell<-sprintf("%s/numcell.dat",exp.folder)
	file_clockwise_eddy<-sprintf("%s/clockwise_eddy.dat",exp.folder)
	file_amplitude_edd<-sprintf("%s/amplitude_eddy.dat",exp.folder)
	file_depth_edd<-sprintf("%s/depth_eddy.dat",exp.folder)
	file_upwelling_edd<-sprintf("%s/upwelling_eddy.dat",exp.folder)
	###
	###


	for (temp_lim in lifetime)
	{
	    print(temp_lim)

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

	#

		#### READING GRID, MASK

	#layer by layer starting from the top
	print(filemask)
	ncm<-nc_open(filemask)
	#varu_3D<-ncm$var[[u_varname]]
	varu_3D<-ncm$var[["lsm"]]
	varu3Dsize<-varu_3D$varsize
	varu3Ddims<-varu_3D$ndims
	print(varu3Ddims)
	print("toto")
	arr_st_3D<-rep(1,varu3Ddims)
	arr_count_3D<-varu3Dsize
	#arr_count_3D[[3]]<-1

	masksuv<-(ncvar_get(ncm,"lsm", start=arr_st_3D, count=arr_count_3D ))>0
	#masksuv<-!is.na(ncvar_get(ncm,u_varname, start=arr_st_3D, count=arr_count_3D ))
	#nc_close(ncm)
	masksuv<- !masksuv   #attension mask is used with 1 associated to land
	gc_clear()


	#filedepth="filedepth.nc"
	#ncmd<-nc_open(filedepth)
	#depth_value<-ncvar_get( ncmd,"depth" ) #1 #ncvar_get( ncm,"depth" )
	#nc_close(ncmd)


        ### READING DX
nav_lon<-ncvar_get(ncm,"lon")
nav_lat<-ncvar_get(ncm,"lat")
nc_close(ncm)

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


#rm(matrix_lon,matrix_lat,matrix_lon_up,matrix_lat_left,boardcond)
gc_clear()

lat_range<-range(nav_lat)

gc_clear()

dz<-1
depth_value <- 1

dx[masksuv]<-NA
dy[masksuv]<-NA
#estentinon for boundary
dx<-array(rbind(dx[(x_window-ext_x_window+1):x_window,],dx,dx[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))
dy<-array(rbind(dy[(x_window-ext_x_window+1):x_window,],dy,dy[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))
matrix_lon_ext<-array(rbind(matrix_lon[(x_window-ext_x_window+1):x_window,],matrix_lon,matrix_lon[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))
matrix_lat_ext<-array(rbind(matrix_lat[(x_window-ext_x_window+1):x_window,],matrix_lat,matrix_lat[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))

####################################
## CREATING MASK FILES
dimX <- ncdim_def("x", "", 1:dim(nav_lon)[1],create_dimvar=FALSE)
dimY <- ncdim_def("y", "", 1:dim(nav_lat)[2],create_dimvar=FALSE)
#dimZ <- ncdim_def("z", "", depth_value[1:z_window], create_dimvar=FALSE)
#dimZ <- ncdim_def("z", "", 1:depth_value, create_dimvar=FALSE)
#dimT <- ncdim_def("time", "",1:n_times_rel,unlim=TRUE, create_dimvar=FALSE)
dimT <- ncdim_def("time", "",1:n_times_rel, create_dimvar=FALSE)
MASK_ind<- ncvar_def("MASK", "", list(dimX,dimY,dimT), -1, longname="  ", prec="single", compression=4,verbose=FALSE)
pos_x_ind<- ncvar_def("lon_eddy", "", list(dimX,dimY,dimT), -1, longname="  ", prec="single", compression=4)
pos_y_ind<- ncvar_def("lat_eddy", "", list(dimX,dimY,dimT), -1, longname="  ", prec="single", compression=4)
radius_ind<- ncvar_def("radius_eddy", "", list(dimX,dimY,dimT), -1, longname="  ", prec="single", compression=4)

lon_ind<- ncvar_def("nav_lon", "", list(dimX,dimY), -1, longname="  ", prec="single", compression=4)
lat_ind<- ncvar_def("nav_lat", "", list(dimX,dimY), -1, longname="  ", prec="single", compression=4)

time_ind<- ncvar_def("time_eddy", "", list(dimT), -1, longname="  ", prec="single", compression=4)

#file_out_mas<-sprintf("%s/MASK_lifetime_%d_%s_%s_totframes_%d_%d_CONE.nc",exp.folder,temp_lim,string_startdate,string_enddate,n_times,z_window)
file_out_mas<-sprintf("%s/POS_MASK_lifetime_%d_%s_%s_totframes_%d_%d_CONE.nc",exp.folder,temp_lim,string_startdate,string_enddate,n_times,z_window)
ncout_mas<-nc_create( file_out_mas, list(MASK_ind,pos_x_ind,pos_y_ind,radius_ind,time_ind,lon_ind,lat_ind) )

# Add Lon and Lat
ncvar_put( ncout_mas, lon_ind, nav_lon, verbose=FALSE )
ncvar_put( ncout_mas, lat_ind, nav_lat, verbose=FALSE )

# Add Time
ncvar_put( ncout_mas, time_ind, seqDATE, verbose=FALSE )
ncatt_put( ncout_mas, time_ind, "origin", sprintf("days since 1970-01-01"), prec="double" )

#
#

##
###################################

###################################
#
# reading data
rel_time<-0
rel_time<-rel_time+framread_first

time_saferange<- c(rel_time+1-(temp_lim+1),rel_time+n_times_rel+(temp_lim+1))
time_saferange[c(time_saferange[1]< 1 ,FALSE)]<- 1
time_saferange[c(FALSE ,time_saferange[2]>n_times)]<- n_times
time_safe_interval<-time_saferange[2]-time_saferange[1]+1


eddies_tot_num<-array(NA,dim=c(n_times))
row_data_index_ed<-scan(file_eddy_index)
row_data_Ekin_eddy_rel_id<-scan(file_Ekin_edd_rel_id)
row_data_median_x<-scan(file_median_x)
row_data_median_y<-scan(file_median_y)
row_data_clock<-scan(file_clockwise_eddy)
row_data_depth<-scan(file_depth_edd)


#--- algorithm for reading files
data_n<-0
for (i in 1:n_times)
{
    data_n<-data_n+1
    eddies_tot_num[i]<-row_data_index_ed[data_n]
    data_n<-data_n+eddies_tot_num[i]
}

max_eddies<-max(eddies_tot_num)
index_ed_t<-array(NA,dim=c(n_times,max_eddies))
Ekin_eddy_rel_id_t<-array(NA,dim=c(n_times,max_eddies))
rel_time_t<-array(NA,dim=c(n_times,max_eddies))
median_x_t<-array(NA,dim=c(n_times,max_eddies))
median_y_t<-array(NA,dim=c(n_times,max_eddies))
clock_t<-array(NA,dim=c(n_times,max_eddies))
depth_t<-array(NA,dim=c(n_times,max_eddies))

data_n<-0
i_rel<-0
for (i in 1:n_times)
{
    data_n<-data_n+1
    if (i >= time_saferange[1] & i <= (time_saferange[2]) & eddies_tot_num[i]>0 ){
        i_rel<-i_rel+1
        index_ed_t[i,1:eddies_tot_num[i]]<-row_data_index_ed[(data_n+1):(data_n+eddies_tot_num[i])]
        Ekin_eddy_rel_id_t[i,1:eddies_tot_num[i]]<-row_data_Ekin_eddy_rel_id[(data_n+1):(data_n+eddies_tot_num[i])]
        median_x_t[i,1:eddies_tot_num[i]]<-row_data_median_x[(data_n+1):(data_n+eddies_tot_num[i])]
        median_y_t[i,1:eddies_tot_num[i]]<-row_data_median_y[(data_n+1):(data_n+eddies_tot_num[i])]
        clock_t[i,1:eddies_tot_num[i]]<-row_data_clock[(data_n+1):(data_n+eddies_tot_num[i])]
        depth_t[i,1:eddies_tot_num[i]]<-row_data_depth[(data_n+1):(data_n+eddies_tot_num[i])]
        rel_time_t[i,1:eddies_tot_num[i]]<-i
    }

    data_n<-data_n+eddies_tot_num[i]
}
##
##to restore real position
median_x_t<-median_x_t + ext_x_window

gc_clear()

#------------- ALL INPUT FILES ARE READ

#------------ ANALYSIS OF THE RESULTS

# grouping values for each eddy

id<-array(index_ed_t[!is.na(index_ed_t)])
#rm(index_ed_t)
Ekin_eddy_rel_id<-array(Ekin_eddy_rel_id_t[!is.na(Ekin_eddy_rel_id_t)])
rm(Ekin_eddy_rel_id_t)
pos_x<-array(median_x_t[!is.na(median_x_t)])
rm(median_x_t)
pos_y<-array(median_y_t[!is.na(median_y_t)])
rm(median_y_t)
clock<-array(clock_t[!is.na(clock_t)])
rm(clock_t)
depth_arr<-array(depth_t[!is.na(depth_t)])
rm(depth_t)
rel_time_arr<-array(rel_time_t[!is.na(rel_time_t)])
rm(rel_time_t)

idf<-factor(id)

eddies_id<-as.integer(levels(idf))
Ekin_eddy_rel_id_list<-array(tapply(Ekin_eddy_rel_id,idf,c))
rel_time_list<-array(tapply(rel_time_arr,idf,c))
pos_x_list<-array(tapply(pos_x,idf,c))
rm(pos_x,rel_time_arr)
pos_y_list<-array(tapply(pos_y,idf,c))
rm(pos_y)
depth_list<-array(tapply(depth_arr,idf,c))
rm(depth_arr)
clock_list<-array(tapply(clock,idf,c))
rm(clock)

gc_clear()

Ekin_eddy_rel_id_array<-array(NA,dim=c(time_safe_interval,length(eddies_id)))
rel_time_array<-array(NA,dim=c(time_safe_interval,length(eddies_id)))
step2_time_array<-array(NA,dim=c(time_safe_interval,length(eddies_id)))
absolute_id_array<-array(NA,dim=c(time_safe_interval,length(eddies_id)))
pos_x_array<-array(NA,dim=c(time_safe_interval,length(eddies_id)))
pos_y_array<-array(NA,dim=c(time_safe_interval,length(eddies_id)))
clock_array<-array(NA,dim=c(time_safe_interval,length(eddies_id)))
depth_array<-array(NA,dim=c(time_safe_interval,length(eddies_id)))
#n_eddies_per_time_length<-array(0,dim=c(n_times))


eddy_vel_x_array<-array(NA,dim=c(time_safe_interval,length(eddies_id)))
eddy_vel_y_array<-array(NA,dim=c(time_safe_interval,length(eddies_id)))
eddy_velmean_x_array<-array(NA,dim=c(time_safe_interval,length(eddies_id)))
eddy_velmean_y_array<-array(NA,dim=c(time_safe_interval,length(eddies_id)))


for (edd in 1:length(eddies_id)  )
{
    temp_Ekin_eddy_rel_id<-array(unlist(Ekin_eddy_rel_id_list[edd]))
    temp_rel_time<-array(unlist(rel_time_list[edd]))
    temp_x<-array(unlist(pos_x_list[edd]))
    temp_y<-array(unlist(pos_y_list[edd]))
    temp_clock<-array(unlist(clock_list[edd]))
    temp_depth<-array(unlist(depth_list[edd]))

    leng<-length(temp_Ekin_eddy_rel_id)

    time_interval<-max(temp_rel_time)-min(temp_rel_time)+1
    if (time_interval<temp_lim) next

    ord<-order(temp_rel_time)
    Ekin_eddy_rel_id_array[1:leng,edd]<-temp_Ekin_eddy_rel_id[ord]
    rel_time_array[1:leng,edd]<-temp_rel_time[ord]
    pos_x_array[1:leng,edd]<-temp_x[ord]
    pos_y_array[1:leng,edd]<-temp_y[ord]
    clock_array[1:leng,edd]<-temp_clock[ord]
    depth_array[1:leng,edd]<-temp_depth[ord]

    temp_l1<-c(array(rel_time_array[2:leng,edd]-rel_time_array[1:(leng-1),edd]),1)
    temp_p1<-c(1,array(rel_time_array[2:leng,edd]-rel_time_array[1:(leng-1),edd]))

    step2_time_array[1:leng,edd]<-(temp_l1+temp_p1)/2

    pos_edd<-cbind(floor(pos_x_array[1:(leng-1),edd]),floor(pos_y_array[1:(leng-1),edd]))
    pos_edd[pos_edd==0]<-1
    eddy_vel_x_array[1:(leng-1),edd]<-(pos_x_array[2:leng,edd]-pos_x_array[1:(leng-1),edd])*dx[pos_edd]/(step2_time_array[1:(leng-1),edd]*days_each_record*24*36)#abs #cm meter

    eddy_vel_y_array[1:(leng-1),edd]<-(pos_y_array[2:leng,edd]-pos_y_array[1:(leng-1),edd])*dy[pos_edd]/(step2_time_array[1:(leng-1),edd]*days_each_record*24*36)#abs #cm meter

    eddy_vel_x_array[leng,edd]<-eddy_vel_x_array[leng-1,edd]
    eddy_vel_y_array[leng,edd]<-eddy_vel_y_array[leng-1,edd]
    eddy_velmean_x_array[1:leng,edd]<-mean(eddy_vel_x_array[1:(leng-1),edd]) #cm/s
    eddy_velmean_y_array[1:leng,edd]<-mean(eddy_vel_y_array[1:(leng-1),edd]) #cm/s

    absolute_id_array[1:leng,edd]<-eddies_id[edd]
}
# gathering values for each eddy

# eliminating row not analysed (lifetime less than the limit
Ekin_eddy_rel_id_array<- Ekin_eddy_rel_id_array[,!is.na(Ekin_eddy_rel_id_array[1,])]
rel_time_array<-rel_time_array[,!is.na(rel_time_array[1,])]
pos_x_array<- pos_x_array[,!is.na(pos_x_array[1,])]
pos_y_array<- pos_y_array[,!is.na(pos_y_array[1,])]
clock_array<-clock_array[,!is.na(clock_array[1,])]
depth_array<-depth_array[,!is.na(depth_array[1,])]
step2_time_array<-step2_time_array[,!is.na(step2_time_array[1,])]
absolute_id_array<- absolute_id_array[,!is.na(absolute_id_array[1,])]
eddy_vel_y_array<-eddy_vel_y_array[,!is.na(eddy_vel_y_array[1,])]
eddy_vel_x_array<-eddy_vel_x_array[,!is.na(eddy_vel_x_array[1,])]
eddy_velmean_y_array<-eddy_velmean_y_array[,!is.na(eddy_velmean_y_array[1,])]
eddy_velmean_x_array<-eddy_velmean_x_array[,!is.na(eddy_velmean_x_array[1,])]

id<-array(rel_time_array[!is.na(rel_time_array)])
Ekin_eddy_rel_id_array<-Ekin_eddy_rel_id_array[!is.na(Ekin_eddy_rel_id_array)]
step2_time_array<-step2_time_array[!is.na(step2_time_array)]
pos_x_array<-pos_x_array[!is.na(pos_x_array)]
pos_y_array<-pos_y_array[!is.na(pos_y_array)]
eddy_vel_x_array<-eddy_vel_x_array[!is.na(eddy_vel_x_array)]
eddy_vel_y_array<-eddy_vel_y_array[!is.na(eddy_vel_y_array)]
eddy_velmean_y_array<-eddy_velmean_y_array[!is.na(eddy_velmean_y_array)]
eddy_velmean_x_array<-eddy_velmean_x_array[!is.na(eddy_velmean_x_array)]
clock_array<-clock_array[!is.na(clock_array)]
depth_array<-depth_array[!is.na(depth_array)]
absolute_id_array<-absolute_id_array[!is.na(absolute_id_array)]

neddies<-length(unique(array(absolute_id_array[!is.na(absolute_id_array)])))


time_focused<- (rel_time):(rel_time+n_times_rel-1)
id_mask<- id %in% time_focused
idf<-factor(id[id_mask])
time_tocount<-as.integer(array(unique(idf)))


rel_index_each_t_list<-array(tapply(array(Ekin_eddy_rel_id_array[id_mask]),idf,c))
step2_time_array_t_list<-array(tapply(array(step2_time_array[id_mask]),idf,c))

pos_x_array_t_list<-array(tapply(array(pos_x_array[id_mask]),idf,c))
pos_y_array_t_list<-array(tapply(array(pos_y_array[id_mask]),idf,c))

depth_array_t_list<-array(tapply(array(depth_array[id_mask]),idf,c))

eddy_vel_x_array_t_list<-array(tapply(array(eddy_vel_x_array[id_mask]),idf,c))
eddy_vel_y_array_t_list<-array(tapply(array(eddy_vel_y_array[id_mask]),idf,c))

eddy_velmean_x_array_t_list<-array(tapply(array(eddy_velmean_x_array[id_mask]),idf,c))
eddy_velmean_y_array_t_list<-array(tapply(array(eddy_velmean_y_array[id_mask]),idf,c))

clock_array_t_list<-array(tapply(array(clock_array[id_mask]),idf,c))
absolute_index_each_t_list<-array(tapply(array(absolute_id_array[id_mask]),idf,c))


rel_index_each_time<-array(NA,dim=c(length(time_tocount),neddies))
step2_time_each_time<-array(NA,dim=c(length(time_tocount),neddies))
pos_x_each_time<-array(NA,dim=c(length(time_tocount),neddies))
pos_y_each_time<-array(NA,dim=c(length(time_tocount),neddies))
eddy_vel_x_each_time<-array(NA,dim=c(length(time_tocount),neddies))
eddy_vel_y_each_time<-array(NA,dim=c(length(time_tocount),neddies))
eddy_velmean_x_each_time<-array(NA,dim=c(length(time_tocount),neddies))
eddy_velmean_y_each_time<-array(NA,dim=c(length(time_tocount),neddies))
clock_each_time<-array(NA,dim=c(length(time_tocount),neddies))
absolute_index_each_time<-array(NA,dim=c(length(time_tocount),neddies))
relative_time_each_time<-array(NA,dim=c(length(time_tocount),neddies))
depth_each_time<-array(NA,dim=c(length(time_tocount),neddies))

if(length(time_tocount)>0){
    for (rel_time in 1:length(time_tocount))
    {
        temp<-array(unlist(rel_index_each_t_list[rel_time]))
        len<-length(temp)
        rel_index_each_time[rel_time,1:len]<-array(unlist(rel_index_each_t_list[rel_time]))
        step2_time_each_time[rel_time,1:len]<-array(unlist(step2_time_array_t_list[rel_time]))
        pos_x_each_time[rel_time,1:len]<-array(unlist(pos_x_array_t_list[rel_time]))
        pos_y_each_time[rel_time,1:len]<-array(unlist(pos_y_array_t_list[rel_time]))
        eddy_vel_x_each_time[rel_time,1:len]<-array(unlist(eddy_vel_x_array_t_list[rel_time]))
        eddy_vel_y_each_time[rel_time,1:len]<-array(unlist(eddy_vel_y_array_t_list[rel_time]))
        eddy_velmean_x_each_time[rel_time,1:len]<-array(unlist(eddy_velmean_x_array_t_list[rel_time]))
        eddy_velmean_y_each_time[rel_time,1:len]<-array(unlist(eddy_velmean_y_array_t_list[rel_time]))
        clock_each_time[rel_time,1:len]<-array(unlist(clock_array_t_list[rel_time]))
        depth_each_time[rel_time,1:len]<-array(unlist(depth_array_t_list[rel_time]))
        absolute_index_each_time[rel_time,1:len]<-array(unlist(absolute_index_each_t_list[rel_time]))
        relative_time_each_time[rel_time,1:len]<-time_tocount[rel_time]
    }
}

rm(id,idf,rel_index_each_t_list,temp,Ekin_eddy_rel_id,Ekin_eddy_rel_id_list,ord,eddies_tot_num,absolute_id_array)
rm(rel_time_list,rel_time_array,temp_Ekin_eddy_rel_id,temp_rel_time,eddies_id,clock_array_t_list,absolute_index_each_t_list)
gc_clear()
rm(step2_time_array_t_list,eddy_vel_x_array_t_list,eddy_vel_y_array_t_list,time_tocount,eddy_velmean_x_array_t_list,eddy_velmean_y_array_t_list,depth_array_t_list)
gc_clear()
#starting
rm(Ekin_eddy_rel_id_array,step2_time_array,eddy_vel_x_array,eddy_vel_y_array,clock_array,eddy_velmean_x_array,eddy_velmean_y_array,depth_array)
gc_clear()

#
#
#################################
#
#STORING RESULTS
print_out<-sprintf("mergeing detection and tracking ..")
print(print_out)


Ekin_rel_surf<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))
Ekin_abs_surf<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))
SLA<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))
Wtot_sum_surf<-array(0,c(ext_x_window+x_window+ext_x_window,y_window))
Wtot_mult_sum<-array(0,c(ext_x_window+x_window+ext_x_window,y_window))
clock_surf<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))

datat_anomal_cy_sum<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))
datat_anomal_anticy_sum<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))
datat_noeddy_sum<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))
mht_anom_cy_sum<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))
zht_anom_cy_sum<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))
mht_anom_anticy_sum<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))
zht_anom_anticy_sum<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))

vy_geo_sum<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))
vy_cy_traj_sum<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))
vy_anticy_traj_sum<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))
vx_geo_sum<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))
vx_cy_traj_sum<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))
vx_anticy_traj_sum<-array(0.,c(ext_x_window+x_window+ext_x_window,y_window))

# Array for the computation of rotational velocity

XINDEX_matr<-array(1:(ext_x_window+x_window+ext_x_window),dim=c(ext_x_window+x_window+ext_x_window,y_window))
YINDEX_matr<-t(array(1:(y_window),dim=c(y_window,ext_x_window+x_window+ext_x_window)))
x_med_matrix<-array(0,dim=c(ext_x_window+x_window+ext_x_window,y_window))
y_med_matrix<-array(0,dim=c(ext_x_window+x_window+ext_x_window,y_window))
sum_ROTVEL_matrix_surf<-array(0,dim=c(ext_x_window+x_window+ext_x_window,y_window))
temp_sum_ROTVEL_matrix<-array(0,dim=c(ext_x_window+x_window+ext_x_window,y_window))
num_cycl_anticycl<-array(0,dim=c(ext_x_window+x_window+ext_x_window,y_window))
dataT_matrix_surf<-array(0,dim=c(ext_x_window+x_window+ext_x_window,y_window))
dataT_meanTpercentage_matrix_surf<-array(0,dim=c(ext_x_window+x_window+ext_x_window,y_window))
vgeo_matrix<-array(0,dim=c(ext_x_window+x_window+ext_x_window,y_window))
ugeo_matrix<-array(0,dim=c(ext_x_window+x_window+ext_x_window,y_window))
#
# relative frame index and eddy index
max_id_prev_edd<-0
rel_time<-0
rel_time<-rel_time+framread_first-1
max_id_prev_edd<-0

TIME_LOCAL<-0


    fileUMEAN<-sprintf("%s/%s_MEAN.nc",dirMEAN,ex)
    fileVMEAN<-sprintf("%s/%s_MEAN.nc",dirMEAN,ex)

#    fileTMEAN<-sprintf("%s/NEATL36_1d_gridTMEAN.nc",dirMEAN)


ncum<-nc_open(fileUMEAN)
ncvm<-nc_open(fileVMEAN)
#nctm<-nc_open(fileTMEAN)

current_date<-startdate
#step_time <- 1 # days

while(current_date<=enddate)
{

count_record<- as.numeric(difftime(current_date,startdate,"days"))+1

string_date<-gsub("-","",current_date)
#next iteration
current_date<- as.Date(current_date + step_time,format="%Y%m%d") #next iteration

print_out<-sprintf("date %s",string_date)
print(print_out)

command<-sprintf("ls -1 %s/%s_%s.nc",dirONE,string_date,ex_lab)
fileSSH<-system(command,inter=TRUE)
print(fileSSH)
fileU<-fileSSH
fileV<-fileU
#NEW AC
command<-sprintf("ls -1 %s/%s_%s.nc",dirSST,string_date,"SST_DYN")

check_exp<-length(grep("VARDYN",ex))>0

if (check_exp) {
 command<-sprintf("ls -1 %s/%s_%s.nc",dirSST,string_date,"SSH_DYN")
}


fileT<-system(command,inter=TRUE)
#fileT<-fileU
nct<-nc_open(fileT)

print_out<-sprintf("reading mask date %s",string_date)
print(print_out)
command<-sprintf("ls -1 %s/Detected_eddies_%s.nc",exp.folder,string_date)
fileMASK<-system(command,inter=TRUE)
print(fileMASK)

#fileS<-fileT

    ncs<-nc_open(fileSSH)
    vars<-ncs$var[[ssh_varname]]
    varsize<-vars$varsize
    vardims<-vars$ndims
    arr_st<-rep(1,vardims)
    arr_count<-varsize
 #   arr_count[[vardims]]<- 1  #one by one record


    ncu<-nc_open(fileU)
    ncv<-nc_open(fileV)
    nct<-nc_open(fileT)
    ncMASK<-nc_open(fileMASK)

    #ncsal<-nc_open(fileS)
#    ncMLD<-nc_open(fileMLD)

    #varu_3D<-ncu$var[["ugosa"]]
    varu_3D<-ncu$var[[u_varname]]
    varu3Dsize<-varu_3D$varsize
    varu3Ddims<-varu_3D$ndims
    arr_st_3D<-rep(1,varu3Ddims)
    arr_count_3D<-varu3Dsize
    #arr_count_3D[[4]]<-1 #one by one record
    #arr_count_3D[[3]]<-1 #one by one level

rel_time<-rel_time+1

start.time<-Sys.time()

#setting the next record, re-initialize the levels
arr_st[vardims]<-1
arr_st_3D[[varu3Ddims]]<-1
#arr_st_3D[[3]]<-1

#setting the mask index
arr_st_3D_dz<- arr_st_3D
#arr_st_3D_dz[[4]]<-1

datas<-ncvar_get( ncs,ssh_varname, start=arr_st, count=arr_count ) #ssh(i,j)
print("ok")
datau<-ncvar_get( ncu,u_varname, start=arr_st_3D, count=arr_count_3D ) #u(i+1/2,j)
datav<-ncvar_get( ncv,v_varname, start=arr_st_3D, count=arr_count_3D ) #v(i,j+1/2)
print("ok again")
#NEW AC
##datat<-ncvar_get( nct,t_varname, start=arr_st_3D, count=arr_count_3D ) #v(i,j+1/2)

##Wtot<-ncvar_get( ncMASK,"Wtot", start=arr_st_3D, count=arr_count_3D ) #v(i,j+1/2)
datat<-ncvar_get( nct,t_varname ) #v(i,j+1/2)
Wtot<-ncvar_get( ncMASK,"Wtot" ) #v(i,j+1/2)
#Wtot<-array(as.integer(Wtot),dim=dim(datat))
print("ok 2again")

#datas<-ncvar_get( ncs,"sla", start=arr_st, count=arr_count ) #ssh(i,j)
#print("ok")
#datau<-ncvar_get( ncu,"ugosa", start=arr_st_3D, count=arr_count_3D ) #u(i+1/2,j)
#datav<-ncvar_get( ncv,"vgosa", start=arr_st_3D, count=arr_count_3D ) #v(i,j+1/2)

#dataMLD<-ncvar_get( ncMLD,"mlotst", start=arr_st, count=arr_count ) #ssh(i,j)
#masksuv<- !is.na(ncvar_get(ncm,"uo", start=arr_st_3D_dz, count=arr_count_3D ))
#NEW AC
#dataMLD<-array(lim_lev,dim(datav))
masksuv<-!is.na(datau)
masksuv<- !masksuv

datas[masksuv]<-NA
datau[masksuv]<-NA
datav[masksuv]<-NA
datat[masksuv]<-NA
#NEW AC
#dataMLD[masksuv]<-NA

#boundary
datas<-array(rbind(datas[(x_window-ext_x_window+1):x_window,],datas,datas[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))
datau<-array(rbind(datau[(x_window-ext_x_window+1):x_window,],datau,datau[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))
datav<-array(rbind(datav[(x_window-ext_x_window+1):x_window,],datav,datav[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))
#NEW AC
datat<-array(rbind(datat[(x_window-ext_x_window+1):x_window,],datat,datat[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))
Wtot<-array(rbind(Wtot[(x_window-ext_x_window+1):x_window,],Wtot,Wtot[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))
#NEW AC
datas[1,]<-NA;
datas[x_window+2*ext_x_window,]<-NA
datau[1,]<-NA;
datau[x_window+2*ext_x_window,]<-NA
datav[1,]<-NA;
datav[x_window+2*ext_x_window,]<-NA
datat[1,]<-NA;
datat[x_window+2*ext_x_window,]<-NA

Wtot[1,]<-0;
Wtot[x_window+2*ext_x_window,]<-0


#dataMLD<-array(rbind(dataMLD[(x_window-ext_x_window+1):x_window,],dataMLD,dataMLD[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))
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

###NEW AC
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

Wtot[is.na(Wtot)]<-0
cc<-ConnCompLabel(Wtot)


# for speeding up, 0's are removed
good_cc<-!is.na(cc) & cc>0
id<-cc[good_cc]
idf<-factor(array(id))
index_ed<-as.integer(levels(idf)) # ID of each eddy
# Added 2nd May 2020
depth_eddy<-array(1,length(index_ed))

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


arr_st[vardims]<-1
arr_st_3D[[varu3Ddims]]<-1
#arr_st_3D[[3]]<-1
arr_st_3D_dz<- arr_st_3D
#arr_st_3D_dz[[4]]<-1
# Toto
arr_count_3D_dz<-arr_count_3D
#arr_count_3D_dz[[4]]<-1


#datat_anomal<-ncvar_get( ncMASK,"T_anom", start=arr_st_3D, count=arr_count_3D ) #v(i,j+1/2)
#datat_noeddy<-ncvar_get( ncMASK,"temp_noeddy", start=arr_st_3D, count=arr_count_3D ) #v(i,j+1/2)
datat_anomal<-ncvar_get( ncMASK,"T_anom")
datat_noeddy<-ncvar_get( ncMASK,"temp_noeddy")
print("ok 3again")

#dataum<-ncvar_get( ncum,u_varname, start=arr_st_3D_dz, count=arr_count_3D_dz ) #u(i+1/2,j)
#datavm<-ncvar_get( ncvm,v_varname, start=arr_st_3D_dz, count=arr_count_3D_dz ) #v(i,j+1/2)
dataum<-ncvar_get( ncum,u_varname ) #u(i+1/2,j)
datavm<-ncvar_get( ncvm,v_varname ) #v(i,j+1/2)
print("ok 4again")

datat_anomal[masksuv]<-NA
dataum[masksuv]<-NA
datavm[masksuv]<-NA
datat_noeddy[masksuv]<-NA
dataum<-array(rbind(dataum[(x_window-ext_x_window+1):x_window,],dataum,dataum[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))
datavm<-array(rbind(datavm[(x_window-ext_x_window+1):x_window,],datavm,datavm[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))
datat_anomal<-array(rbind(datat_anomal[(x_window-ext_x_window+1):x_window,],datat_anomal,datat_anomal[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))
datat_noeddy<-array(rbind(datat_noeddy[(x_window-ext_x_window+1):x_window,],datat_noeddy,datat_noeddy[1:ext_x_window,]),dim=c(x_window+2*ext_x_window,y_window))

dataum[1,]<-NA;
dataum[x_window+2*ext_x_window,]<-NA
datavm[1,]<-NA;
datavm[x_window+2*ext_x_window,]<-NA

datat_noeddy[1,]<-NA;
datat_noeddy[x_window+2*ext_x_window,]<-NA
datat_anomal[1,]<-NA;
datat_anomal[x_window+2*ext_x_window,]<-NA
 
### calculating no-eddy features

temp<- (    (datau-dataum)*(datau-dataum)+(datav-datavm)*(datav-datavm)   ) /2.
temp[is.na(temp)]<-0
Ekin_abs_surf<-Ekin_abs_surf+temp



index_ed_t_arr1<-index_ed_t[rel_time,]


if( sum( rel_time %in% relative_time_each_time ,na.rm=TRUE) < 1 ){

    datat_noeddy_sum<-datat_noeddy_sum+datat

    end.time<-Sys.time()
    exec.time<-end.time-start.time
    cat("no eddy found rel time=",rel_time)
    print(exec.time)

    ncvar_put( ncout_mas, MASK_ind, Wtot[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=c(1,1, count_record), count=c(x_window,y_window,1), verbose=FALSE )

    next
}

TIME_LOCAL<-TIME_LOCAL+1
temp_x<-pos_x_each_time[TIME_LOCAL,]
temp_y<-pos_y_each_time[TIME_LOCAL,]
temp_rel<-rel_index_each_time[TIME_LOCAL,]
temp_step2<-step2_time_each_time[TIME_LOCAL,]
temp_x_vel<-eddy_vel_x_each_time[TIME_LOCAL,]
temp_y_vel<-eddy_vel_y_each_time[TIME_LOCAL,]
temp_x_velmean<-eddy_velmean_x_each_time[TIME_LOCAL,]
temp_y_velmean<-eddy_velmean_y_each_time[TIME_LOCAL,]
temp_clock<-clock_each_time[TIME_LOCAL,]
temp_DEPTH<-depth_each_time[TIME_LOCAL,]
temp_absol_id<-absolute_index_each_time[TIME_LOCAL,]

temp_rel<-temp_rel[!is.na(temp_rel)]
temp_x<-temp_x[!is.na(temp_x)]
temp_y<-temp_y[!is.na(temp_y)]
temp_step2<-temp_step2[!is.na(temp_step2)]
temp_x_vel<-temp_x_vel[!is.na(temp_x_vel)]
temp_y_vel<-temp_y_vel[!is.na(temp_y_vel)]
temp_x_velmean<-temp_x_velmean[!is.na(temp_x_velmean)]
temp_y_velmean<-temp_y_velmean[!is.na(temp_y_velmean)]
temp_clock<-temp_clock[!is.na(temp_clock)]
temp_DEPTH<-temp_DEPTH[!is.na(temp_DEPTH)]
temp_absol_id<-temp_absol_id[!is.na(temp_absol_id)]


Wtot[!(cc %in% temp_rel)]<-0
cc[!(cc %in% temp_rel)]<-NA

good_cc<-!is.na(cc) & cc>0

len_lev<-length(levels(factor(array(cc[good_cc]))))
idf<-factor(array(cc[good_cc]),labels=1:(len_lev))


datat_noeddy[!Wtot]<-datat[!Wtot]
datat_noeddy_sum<-datat_noeddy_sum+datat_noeddy

dist_temp<-round(sqrt(x_matr^2+y_matr^2))
#


index_ed_t_arr1[!(index_ed_t_arr1 %in% temp_absol_id)]<-NA
index_ed_t_arr1<-index_ed_t_arr1[!is.na(index_ed_t_arr1)]

NUM_eddy<-index_ed_t_arr1[as.integer(idf)]
NUM_eddy_matr<-array(NA,dim=dim(datas))
NUM_eddy_matr[as.logical(Wtot)]<-NUM_eddy
#######
len_lev_n<-length(levels(factor(array(cc[good_cc]))))
idf_n<-factor(array(cc[good_cc]),labels=1:(len_lev_n))
index_ed<-as.integer(levels(idf_n))#to get all the eddies
xy_Wtot<-which(good_cc,arr.ind=TRUE)
median_x<-array(tapply(xy_Wtot[,1],idf_n,mean))
median_y<-array(tapply(xy_Wtot[,2],idf_n,mean))
#cycl_eddy1<-array(tapply(array(Wmask1[good_cc]),idf_n,sum))#clockwise
#cycl_eddy2<-array(tapply(-array(Wmask2[good_cc]),idf_n,sum))#counterclock
#cycl_eddy<-cycl_eddy1+cycl_eddy2
#cycl_eddy[cycl_eddy>0]<- 1 #do not switch!!
#cycl_eddy[cycl_eddy<0]<- -1
sum_ROTVEL<-array(tapply(array(sum_ROTVEL_matrix_surf[good_cc]),idf_n,mean,na.rm=TRUE))
cycl_eddy<-sum_ROTVEL
cycl_eddy[cycl_eddy>0]<- 1 #do not switch!!
cycl_eddy[cycl_eddy<0]<- -1

######
# most important indeed to be the same as new
ord<-order(temp_rel)

temp_x<-temp_x[ord]
temp_y<-temp_y[ord]
temp_x_vel<-temp_x_vel[ord]
temp_y_vel<-temp_y_vel[ord]
temp_step2<-temp_step2[ord]
temp_clock<-temp_clock[ord]
temp_x_velmean<-temp_x_velmean[ord]
temp_y_velmean<-temp_y_velmean[ord]
temp_DEPTH<-temp_DEPTH[ord]

Wtot<-as.logical(Wtot)
Wtot<-array(Wtot,dim(datau))

#remove boudary
ttmep<-median_x[idf]
mediax_matrix<-array(0,dim=dim(Wtot))
mediax_matrix[Wtot]<-ttmep

Wtot_no_boundary<-Wtot
Wtot_no_boundary[!(mediax_matrix>(ext_x_window) & mediax_matrix<(ext_x_window+x_window))]<- FALSE


cc_noboudary<-ConnCompLabel(Wtot_no_boundary)
good_cc_noboudary<-!is.na(cc_noboudary) & cc_noboudary>0
len_lev_noboudary<-length(levels(factor(array(cc_noboudary[good_cc_noboudary]))))
idf_noboudary<-factor(array(cc_noboudary[good_cc_noboudary]),labels=1:(len_lev_noboudary))

#remove boudary

temp<-temp_clock[idf_noboudary]
clock_matrix<-array(0,dim=dim(Wtot_no_boundary))
clock_matrix[Wtot_no_boundary]<-temp

clock_full_matrix<-array(0,dim=dim(Wtot_no_boundary))
temp<-temp_clock[as.integer(idf_noboudary)]
clock_full_matrix[Wtot_no_boundary]<- temp
clock_full_matrix[!Wtot_no_boundary]<-NA

clock_12<-clock_full_matrix
clock_12[clock_12< 0]<- 2

datat_anomal_cy_sum<-datat_anomal_cy_sum+datat_anomal*Wtot*(clock_full_matrix>0 & !is.na(clock_full_matrix))
datat_anomal_anticy_sum<-datat_anomal_anticy_sum+datat_anomal*Wtot*(clock_full_matrix<0 & !is.na(clock_full_matrix))

temp<-temp_x_vel[idf_noboudary]
xvel_matrix<-array(0,dim=dim(Wtot_no_boundary))
xvel_matrix[Wtot_no_boundary]<-temp
temp<-temp_y_vel[idf_noboudary]
yvel_matrix<-array(0,dim=dim(Wtot_no_boundary))
yvel_matrix[Wtot_no_boundary]<-temp

temp<-temp_x_velmean[idf_noboudary]
xvelmean_matrix<-array(0,dim=dim(Wtot_no_boundary))
xvelmean_matrix[Wtot_no_boundary]<-temp
temp<-temp_y_velmean[idf_noboudary]
yvelmean_matrix<-array(0,dim=dim(Wtot_no_boundary))
yvelmean_matrix[Wtot_no_boundary]<-temp


mht_anom_cy_sum<-mht_anom_cy_sum+datat_anomal*Wtot* yvelmean_matrix*(clock_full_matrix>0 & !is.na(clock_full_matrix))
zht_anom_cy_sum<-zht_anom_cy_sum+datat_anomal*Wtot* xvelmean_matrix*(clock_full_matrix>0 & !is.na(clock_full_matrix))
mht_anom_anticy_sum<-mht_anom_anticy_sum+datat_anomal*Wtot* yvelmean_matrix*(clock_full_matrix<0 & !is.na(clock_full_matrix))
zht_anom_anticy_sum<-zht_anom_anticy_sum+datat_anomal*Wtot* xvelmean_matrix*(clock_full_matrix<0 & !is.na(clock_full_matrix))

#mht_anom_cy_sum<-mht_anom_cy_sum+datat_anomal*Wtot* yvel_matrix*(clock_full_matrix>0 & !is.na(clock_full_matrix))
#zht_anom_cy_sum<-zht_anom_cy_sum+datat_anomal*Wtot* xvel_matrix*(clock_full_matrix>0 & !is.na(clock_full_matrix))
#mht_anom_anticy_sum<-mht_anom_anticy_sum+datat_anomal*Wtot* yvel_matrix*(clock_full_matrix<0 & !is.na(clock_full_matrix))
#zht_anom_anticy_sum<-zht_anom_anticy_sum+datat_anomal*Wtot* xvel_matrix*(clock_full_matrix<0 & !is.na(clock_full_matrix))


vx_cy_traj_sum<-vx_cy_traj_sum+Wtot* xvelmean_matrix*(clock_full_matrix>0 & !is.na(clock_full_matrix))
vx_anticy_traj_sum<-vx_anticy_traj_sum+Wtot* xvelmean_matrix*(clock_full_matrix<0 & !is.na(clock_full_matrix))

vy_cy_traj_sum<-vy_cy_traj_sum+Wtot*yvelmean_matrix*(clock_full_matrix>0 & !is.na(clock_full_matrix))
vy_anticy_traj_sum<-vy_anticy_traj_sum+Wtot*yvelmean_matrix*(clock_full_matrix<0 & !is.na(clock_full_matrix))

#vx_cy_traj_sum<-vx_cy_traj_sum+Wtot* xvel_matrix*(clock_full_matrix>0 & !is.na(clock_full_matrix))
#vx_anticy_traj_sum<-vx_anticy_traj_sum+Wtot* xvel_matrix*(clock_full_matrix<0 & !is.na(clock_full_matrix))

#vy_cy_traj_sum<-vy_cy_traj_sum+Wtot*yvel_matrix*(clock_full_matrix>0 & !is.na(clock_full_matrix))
#vy_anticy_traj_sum<-vy_anticy_traj_sum+Wtot*yvel_matrix*(clock_full_matrix<0 & !is.na(clock_full_matrix))

##geostrophic mean on the eddy
ugeo_matrix[,]<-0
sum_datau<-array(tapply(datau[good_cc_noboudary],idf_noboudary,mean,na.rm=TRUE))
temp<-sum_datau[as.integer(idf_noboudary)]
ugeo_matrix[Wtot_no_boundary]<- temp
vx_geo_sum<-vx_geo_sum+ugeo_matrix
#vx_geo_sum<-vx_geo_sum+datau

vgeo_matrix[,]<-0
sum_datav<-array(tapply(datav[good_cc_noboudary],idf_noboudary,mean,na.rm=TRUE))
temp<-sum_datav[as.integer(idf_noboudary)]
vgeo_matrix[Wtot_no_boundary]<- temp
vy_geo_sum<-vy_geo_sum+vgeo_matrix
#vy_geo_sum<-vy_geo_sum+datav


#temp<-depth_value[temp_DEPTH[idf]]
temp<-1
DEPTH_matrix<-array(0,dim=dim(Wtot_no_boundary))
DEPTH_matrix[Wtot_no_boundary]<-temp


temp<-temp_step2[idf_noboudary]
Wtot_mult<-array(1,dim=dim(Wtot_no_boundary))
Wtot_mult[Wtot_no_boundary]<-temp
Wtot_mult_zero<-array(0,dim=dim(Wtot_no_boundary))
Wtot_mult_zero[Wtot_no_boundary]<-temp

temp<-clock_matrix*Wtot_mult
temp[is.na(temp)]<-0
clock_surf<-clock_surf+temp

Wtot_mult_zero[Wtot_mult_zero==1]<-0
#rm(Wtot_multipliaction)
gc_clear()

temp<-( (datau-dataum)*(datau-dataum)  +(datav-datavm)*(datav-datavm)  )*Wtot_mult/2 #*dx*dy*dz_1l
temp[!Wtot_no_boundary]<-0
temp[is.na(temp)]<-0
Ekin_rel_surf<-Ekin_rel_surf+temp


temp<-datas* Wtot_mult
temp[!Wtot_no_boundary]<-0
temp[is.na(temp)]<-0
SLA<-SLA+temp


Wtot_sum_surf<-Wtot_sum_surf+as.integer(Wtot_no_boundary)*Wtot_mult
Wtot_mult_sum<-Wtot_mult_sum+ as.integer(Wtot_mult==1.5 | Wtot_mult==2)

temp_y[temp_y<1]<-1;temp_y[temp_y>y_window]<-y_window;
pos_mean<-array(cbind(round(temp_x),round(temp_y)),dim=c(length(temp_x),2)) #already summed pos + ext_x

num_eddy<-array(0,dim=c(2*ext_x_window+x_window,y_window))
num_eddy[pos_mean]<-1
num_eddy<-num_eddy*Wtot_mult


eddy_mask<-Wtot_no_boundary
n_depth<-1



temp_sum_ROTVEL_matrix[,]<-0
num_cycl_anticycl[,]<-0

temp_sum_ROTVEL<-sum_ROTVEL
temp_sum_ROTVEL[sum_ROTVEL<0]<-0 #cycl

temp_sum_ROTVEL_matrix[pos_mean]<-temp_sum_ROTVEL
num_cycl_anticycl[abs(temp_sum_ROTVEL_matrix)>0]<-1

num_cycl_anticycl<-num_cycl_anticycl*Wtot_mult
temp_sum_ROTVEL_matrix<-temp_sum_ROTVEL_matrix*Wtot_mult

temp_sum_ROTVEL_matrix[,]<-0
num_cycl_anticycl[,]<-0

temp_sum_ROTVEL<-sum_ROTVEL
temp_sum_ROTVEL[sum_ROTVEL>0]<-0 #anticycl

temp_sum_ROTVEL_matrix[pos_mean]<-temp_sum_ROTVEL
num_cycl_anticycl[abs(temp_sum_ROTVEL_matrix)>0]<-1

num_cycl_anticycl<-num_cycl_anticycl*Wtot_mult
temp_sum_ROTVEL_matrix<-temp_sum_ROTVEL_matrix*Wtot_mult


x_med_matrix[,]<-0
y_med_matrix[,]<-0
sum_ROTVEL_matrix_surf[,]<-0
temp<-temp_x[as.integer(idf_noboudary)]
x_med_matrix[Wtot_no_boundary]<- temp
temp<-temp_y[as.integer(idf_noboudary)]
y_med_matrix[Wtot_no_boundary]<- temp

x_matr<-round(XINDEX_matr-x_med_matrix)
y_matr<-round(YINDEX_matr-y_med_matrix)

#x_matr_exteddy<-x_matr
#y_matr_exteddy<-y_matr
#x_matr_exteddy[Wtot]<-NA
#y_matr_exteddy[Wtot]<-NA

x_matr[!Wtot_no_boundary]<-NA
y_matr[!Wtot_no_boundary]<-NA
ROTVEL<-((datau)*y_matr*dy-x_matr*dx*(datav))/(sqrt( (y_matr*dy)^2+(x_matr*dx)^2))
sum_ROTVEL<-array(tapply(ROTVEL[good_cc_noboudary],idf_noboudary,mean,na.rm=TRUE))
temp<-sum_ROTVEL[as.integer(idf_noboudary)]
sum_ROTVEL_matrix_surf[Wtot_no_boundary]<- temp

#lon and lat
lon_array<-array(tapply(matrix_lon_ext[good_cc_noboudary],idf_noboudary,mean))
lat_array<-array(tapply(matrix_lat_ext[good_cc_noboudary],idf_noboudary,mean))

matrix_lon_mask<-array(0,dim=dim(Wtot_no_boundary))
temp<-lon_array[as.integer(idf_noboudary)]
matrix_lon_mask[Wtot_no_boundary]<- temp

matrix_lat_mask<-array(0,dim=dim(Wtot_no_boundary))
temp<-lat_array[as.integer(idf_noboudary)]
matrix_lat_mask[Wtot_no_boundary]<- temp

### radius

Patch<-PatchStat(cc_noboudary[good_cc_noboudary])
num_cell<-Patch$n.cell

factor_dy_id<-array(tapply(dy[good_cc_noboudary],idf_noboudary,mean,na.rm=TRUE))
factor_dx_id<-array(tapply(dx[good_cc_noboudary],idf_noboudary,mean,na.rm=TRUE))

radius<-sqrt(num_cell* factor_dy_id*factor_dx_id/3.14)
matrix_radius<-array(0,dim=dim(Wtot_no_boundary))
temp<-radius[as.integer(idf_noboudary)]
matrix_radius[Wtot_no_boundary]<- temp

# Toto
matrix_radius<-matrix_radius*clock_full_matrix
###

ROTVEL<-((datau)*y_matr*dy-x_matr*dx*(datav))/(sqrt( (y_matr*dy)^2+(x_matr*dx)^2))*Wtot_mult
NUMPOINT<-Wtot_no_boundary*Wtot_mult

maxx<-max(abs(x_matr),na.rm=TRUE)+1
maxy<-max(abs(y_matr),na.rm=TRUE)+1

dist_fix<- (maxx+x_matr)+1000*(maxy+y_matr)+1000000*clock_12


end.time<-Sys.time()
exec.time<-end.time-start.time
cat("rel time=",rel_time)
print(exec.time)

print(dim(Wtot_no_boundary))


Wtot_put<-Wtot_no_boundary[(ext_x_window+1):(x_window+ext_x_window),1:y_window]

#ncvar_put( ncout_mas, MASK_ind, Wtot[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=c(1,1), count=c(x_window,y_window), verbose=TRUE)
#ncvar_put( ncout_mas, MASK_ind, Wtot[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=c(1,1,count_record), count=c(x_window,y_window,1), verbose=TRUE)
ncvar_put( ncout_mas, MASK_ind, Wtot_put, start=c(1,1,count_record), count=c(x_window,y_window,1), verbose=FALSE)
rm(Wtot_put)

#ncvar_put( ncout_mas, pos_x_ind, matrix_lon_mask[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=c(1,1,count_record), count=c(x_window,y_window,1), verbose=FALSE )
#ncvar_put( ncout_mas, pos_y_ind, matrix_lat_mask[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=c(1,1,count_record), count=c(x_window,y_window,1), verbose=FALSE )
#ncvar_put( ncout_mas, radius_ind, matrix_radius[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=c(1,1,count_record), count=c(x_window,y_window,1), verbose=FALSE )

nc_close(ncs)
nc_close(ncu)
nc_close(ncv)
nc_close(nct)
#nc_close(ncM)
#nc_close(ncm)
nc_close(ncMASK)
#nc_close(ncsal)
#nc_close(ncMLD)
}
nc_close(ncum)
nc_close(ncvm)
#nc_close(nctm)
#nc_close(ncm)

gc_clear()


dimX <- ncdim_def("x", "", 1:x_window,create_dimvar=FALSE)
dimY <- ncdim_def("y", "", 1:y_window,create_dimvar=FALSE)
#
Wtot_sum_surf_ind<- ncvar_def("Wtot_surf", "counts", list(dimX,dimY), 9.96921e+36, longname=" number of passages each point durface ", prec="single",compression=4)
EKE_rel_surf_ind<- ncvar_def("EKE_rel_surf", "energy", list(dimX,dimY), 9.96921e+36, longname="  ", prec="single",compression=4)
EKE_tot_surf_ind<- ncvar_def("EKE_tot_surf", "energy", list(dimX,dimY), 9.96921e+36, longname="  ", prec="single",compression=4)
datat_noeddy_ind<- ncvar_def("SST_NOEDDY", "°C", list(dimX,dimY), 9.96921e+36, longname="  ", prec="single",compression=4)
datat_anomaly_cy_ind<- ncvar_def("SST_ANOMALY_CY", "°C", list(dimX,dimY), 9.96921e+36, longname="  ", prec="single",compression=4)
datat_anomaly_anticy_ind<- ncvar_def("SST_ANOMALY_ANTICY", "°C", list(dimX,dimY), 9.96921e+36, longname="  ", prec="single",compression=4)
mht_anom_cy_ind<- ncvar_def("MHT_ANOM_CY", "°C*cm/s", list(dimX,dimY), 9.96921e+36, longname="  ", prec="single",compression=4)
zht_anom_cy_ind<- ncvar_def("ZHT_ANOM_CY", "°C*cm/s", list(dimX,dimY), 9.96921e+36, longname="  ", prec="single",compression=4)
mht_anom_anticy_ind<- ncvar_def("MHT_ANOM_ANTICY", "°C*cm/s", list(dimX,dimY), 9.96921e+36, longname="  ", prec="single",compression=4)
zht_anom_anticy_ind<- ncvar_def("ZHT_ANOM_ANTICY", "°C*cm/s", list(dimX,dimY), 9.96921e+36, longname="  ", prec="single",compression=4)

vy_geo_ind<- ncvar_def("VY_GEOST", "m/s", list(dimX,dimY), 9.96921e+36, longname="  m/s ", prec="single",compression=4)
vy_cy_traj_ind<- ncvar_def("VY_CY_TRAJEC", "m/s", list(dimX,dimY), 9.96921e+36, longname=" multiplied by 0.01 for m/s  ", prec="single",compression=4)
vy_anticy_traj_ind<- ncvar_def("VY_ANTICY_TRAJEC", "m/s", list(dimX,dimY), 9.96921e+36, longname=" multiplied by 0.01 for m/s  ", prec="single",compression=4)
vx_geo_ind<- ncvar_def("VX_GEOST", "m/s", list(dimX,dimY), 9.96921e+36, longname=" m/s   ", prec="single",compression=4)
vx_cy_traj_ind<- ncvar_def("VX_CY_TRAJEC", "m/s", list(dimX,dimY), 9.96921e+36, longname=" multiplied by 0.01 for m/s  ", prec="single",compression=4)
vx_anticy_traj_ind<- ncvar_def("VX_ANTICY_TRAJEC", "m/s", list(dimX,dimY), 9.96921e+36, longname=" multiplied by 0.01 for m/s  ", prec="single",compression=4)

clock_surf_ind<-ncvar_def("clock_surf", "adim", list(dimX,dimY), 9.96921e+36, longname=" clockwise/counter ", prec="single",compression=4)
SLA_ind<-ncvar_def("SSH", "m", list(dimX,dimY), 9.96921e+36, longname="SSH ", prec="single",compression=4)
Wtot_mult_ind<-ncvar_def("Wtotmult", "count", list(dimX,dimY), 9.96921e+36, longname="countmissed ", prec="single",compression=4)
nav_lon_ind<-ncvar_def("nav_lon", "degrees_east", list(dimX,dimY), 9.96921e+36, longname="", prec="single",compression=4)
nav_lat_ind<-ncvar_def("nav_lat", "degrees_north", list(dimX,dimY), 9.96921e+36, longname=" ", prec="single",compression=4)


file_out<-sprintf("%s/Features_eddies_lifetime_%d_%s_%s_totframes_%d_TEMP_INTERP_MAX_DEPTH_%d_CONE.nc",exp.folder,temp_lim,string_startdate,string_enddate,n_times,z_window)
ncout<-nc_create( file_out, list(Wtot_sum_surf_ind,EKE_rel_surf_ind,
clock_surf_ind,EKE_tot_surf_ind,SLA_ind,Wtot_mult_ind,nav_lat_ind,nav_lon_ind,datat_noeddy_ind,datat_anomaly_anticy_ind,datat_anomaly_cy_ind,mht_anom_cy_ind,zht_anom_cy_ind,mht_anom_anticy_ind,zht_anom_anticy_ind,
vy_geo_ind,vy_cy_traj_ind,vy_anticy_traj_ind,vx_geo_ind,vx_cy_traj_ind,vx_anticy_traj_ind),force_v4=TRUE )

start_arr<-c(1,1)
count_arr<-c(x_window,y_window)
ncvar_put( ncout, nav_lat_ind, nav_lat, start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, nav_lon_ind, nav_lon, start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, datat_noeddy_ind, datat_noeddy_sum[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, datat_anomaly_anticy_ind, datat_anomal_anticy_sum[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, datat_anomaly_cy_ind, datat_anomal_cy_sum[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, EKE_rel_surf_ind, Ekin_rel_surf[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, EKE_tot_surf_ind, Ekin_abs_surf[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, zht_anom_cy_ind, zht_anom_cy_sum[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, mht_anom_cy_ind, mht_anom_cy_sum[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, zht_anom_anticy_ind, zht_anom_anticy_sum[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, mht_anom_anticy_ind, mht_anom_anticy_sum[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )

ncvar_put( ncout, vy_geo_ind, vy_geo_sum[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, vy_cy_traj_ind, 0.01*vy_cy_traj_sum[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, vy_anticy_traj_ind, 0.01*vy_anticy_traj_sum[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, vx_geo_ind, vx_geo_sum[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, vx_cy_traj_ind, 0.01*vx_cy_traj_sum[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, vx_anticy_traj_ind, 0.01*vx_anticy_traj_sum[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )

ncvar_put( ncout, clock_surf_ind, clock_surf[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, SLA_ind, SLA[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, Wtot_sum_surf_ind, Wtot_sum_surf[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )
ncvar_put( ncout, Wtot_mult_ind, Wtot_mult_sum[(ext_x_window+1):(x_window+ext_x_window),1:y_window], start=start_arr, count=count_arr, verbose=FALSE )


nc_close(ncout)
nc_close(ncout_mas)

#warnings()

}
