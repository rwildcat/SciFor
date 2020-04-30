! SciFor tests - io_netcdf
! Library: NetCDF
!
! RSolano, Mar/03/2020
!
! To compile the test:
! 		$ gfortran io_netcdf_scifor-test.f90  -logpf `pkg-config --cflags --libs scifor` 

program netcdf_test

   use scifor
   use ogpf
   
   implicit none
   character(256) :: fname, vname, fillname
   character(10), allocatable :: dnames(:)   ! dim names
   character(10)              :: xdim, ydim  ! xdim and y dim names
   integer, allocatable       ::  vdimsl(:)  ! dims lens
   integer              :: i, ndims
   real, allocatable    :: var2d(:,:), var3d(:,:,:)
   real(8), allocatable :: lat(:), lon(:), glat(:,:), glon(:,:)
   real                  :: fillval
   type(gpf)            :: plt

   ! lon = 180 ; 
   ! lat = 170 ;
   ! time = UNLIMITED ; // (24 currently)
   ! double lon(lon) ;
   ! double lat(lat) ;
   ! float tos(time, lat, lon) ; sst, K, fill=1e20
   fname = 'tos_O1_2001-2002.nc'
   vname = 'tos' 
   xdim = 'lon'
   ydim = 'lat'
   fillname = '_FillValue'
   
   ! main

   print *, 'io_netcdf:'
   print *, '   fname = ', trim(fname)
   print *, '   vname = ', trim(vname)
   
   ! get dims
   call io_nc_getvardims(fname, vname, vdimsl, ndims_=ndims, dnames_=dnames)  

   ! show 
   print *, '   dim names : ', trim(vname), ' (', (trim(dnames(i)) // ' ', i=1, size(dnames)), ')'
   print *, '   dim len   : ', trim(vname), ' (', vdimsl, ')'

   ! read lat, lon
   call io_nc_getvar(fname, trim(ydim), lat)
   call io_nc_getvar(fname, trim(xdim), lon)
   !call io_nc_getvar(fname, trim(ydim), lat, start_=86)
   !call io_nc_getvar(fname, trim(xdim), lon, start_=91)

   ! allocate mem for plotting grids (glat,glon)
   allocate( glat(size(lat), size(lon)) )
   allocate( glon(size(lat), size(lon)) )
   
   ! make grids for plotting (meshgrid(X, Y, RANGEX, RANGEY)
   call meshgrid(glon, glat, lon, lat)

   ! get data
   ! call io_nc_getvar(fname, vname, var2d)
   call io_nc_getvar(fname, vname, var3d)
   !call io_nc_getvar(fname, vname, var3d, count_=[1,85,90], start_=[1,86,91])
   
   ! get att fillvalue
   call io_nc_getatt(fname, vname, fillname, fillval)
   print *, trim(fillname), ' = ', fillval

   ! "mask" nodata, use min val
   ! where(var2d==fillval) var2d = minval(var2d)
   where(var3d==fillval) var3d = minval(var3d)

   print *, 'Plotting...'
   call plt%options('set term aqua')
   call plt%options('set view equal xy')
   call plt%options('set view map')

   ! call plt%surf(1d0*log(var2d))
   ! call plt%surf(glon, glat, 1d0*var2d, lspec='with image')
   call plt%surf(glon, glat, 1d0*var3d(1,:,:), lspec='with image')

   ! call plt%animation_start(1)
   ! do i=1, 23
   !   call plt%surf(glon, glat, 1d0*var3d(i,:,:), lspec='with image')   
   ! end do
   ! call plt%animation_show()
      
end program
