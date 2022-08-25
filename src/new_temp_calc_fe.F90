!**********************************************************************************
! Module calculates the new temperature at new time step using Forward Euler method
!**********************************************************************************
module calc_new_T_fe
use Constants
use, intrinsic :: iso_fortran_env
implicit none
public:: calc_new_T_surf_diff, calc_new_T_deep_diff,calc_new_T_surf_notrns, &
        calc_new_T_deep_notrns
private


contains
!*****************************************************************************************
! Subroutines to calculate new temperature using the forward Euler method (No transport)
!*****************************************************************************************
subroutine calc_new_T_surf_notrns(T,i,j,h,s_notrns,F_c,F_a,T_new)
    real(real64),dimension(:,:,:),intent(in) :: T
    real(real64),dimension(:,:,:),intent(out) :: T_new
    real(real64),dimension(:,:), intent(in) :: F_a, F_c
    integer(int64), intent(in) :: h,i,j
    real(real64), intent(in) :: s_notrns

    T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4) + T(h,i,j)

end subroutine calc_new_T_surf_notrns

subroutine calc_new_T_deep_notrns(T,i,j,h,s_notrns,F_c,T_new)
    real(real64),dimension(:,:,:),intent(in) :: T
    real(real64),dimension(:,:,:),intent(out) :: T_new
    integer(int64), intent(in) :: h,i,j
    real(real64),dimension(:,:), intent(in) :: F_c
    real(real64), intent(in) :: s_notrns

    T_new(h,i,j) = s_notrns*(-F_c(i,j))+T(h,i,j)

end subroutine calc_new_T_deep_notrns

!*****************************************************************************************************************************
! At maximum or minimum latitude, the adjacent latitude point above or bellow, respectively, is located at the same latitude,
! but at the new longitudinal point (from lon to new_lon)
!******************************************************************************************************************************
function calculate_new_lon(lon) result(new_lon)
    integer(int64), intent(in) :: lon
    integer(int64) :: new_lon

    new_lon = lon+N_LONS/2
    if (new_lon > N_LONS) then
        new_lon = new_lon - N_LONS
    end if
end function calculate_new_lon

!*****************************************************************************************
! Subroutines to calculate new temperature using the forward Euler method (with diffusion)
!*****************************************************************************************
subroutine calc_new_T_surf_diff(T,i,j,h,d_theta,d_phi,s_notrns,s_dfsn,F_c,F_a,theta,T_new)
    real(real64),dimension(:,:,:),intent(in) :: T
    integer(int64), intent(in) :: h,i,j
    real(real64), intent(in) :: d_theta,d_phi,s_notrns,s_dfsn,theta
    real(real64),dimension(:,:), intent(in) :: F_a, F_c
    real(real64),dimension(:,:,:),intent(out) :: T_new

    if(i==1) then
        ! T(h,i-1,j) = T(h,i,calculate_new_lon(j))
        if(j==1) then
             ! j-1 = N_LONS
            T_new(h,i,j) =  s_notrns*(F_a(i,j)+F_c(i,j)+SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                        s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2 - &
                        tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
                         T(h,i,j)

        elseif(j==N_LONS) then
            ! j+1 = 1
            T_new(h,i,j) =  s_notrns*(F_a(i,j)+F_c(i,j)+SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                        s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2 - &
                        tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                         T(h,i,j)

        else
            T_new(h,i,j) =  s_notrns*(F_a(i,j)+F_c(i,j)+SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                        s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2 - &
                        tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                         T(h,i,j)


        end if

    elseif(i==N_LATS) then
        !T(h,i+1,j) = T(h,i,calculate_new_lon(j))
        if(j==1) then
            ! j-1 = N_LONS
            T_new(h,i,j) =  s_notrns*(F_a(i,j)+F_c(i,j)+SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                        s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2 - &
                        tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
                         T(h,i,j)
        elseif(j==N_LONS) then
            ! j+1 = 1
            T_new(h,i,j) =  s_notrns*(F_a(i,j)+F_c(i,j)+SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                        s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2 - &
                        tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                         T(h,i,j)
        else

            T_new(h,i,j) =  s_notrns*(F_a(i,j)+F_c(i,j)+SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                        s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2 - &
                        tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                         T(h,i,j)
        end if

    else
        if(j==1) then
            ! j-1 = N_LONS
            T_new(h,i,j) =  s_notrns*(F_a(i,j)+F_c(i,j)+SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                        s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2 - &
                        tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
                         T(h,i,j)
        elseif(j==N_LONS) then
            ! j+1 = 1
            T_new(h,i,j) =  s_notrns*(F_a(i,j)+F_c(i,j)+SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                        s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2 - &
                        tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                         T(h,i,j)
        else

            T_new(h,i,j) =  s_notrns*(F_a(i,j)+F_c(i,j)+SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                        s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2 - &
                        tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                         T(h,i,j)
        end if

    endif

end subroutine calc_new_T_surf_diff

subroutine calc_new_T_deep_diff(T,i,j,h,d_theta,d_phi,s_notrns,s_dfsn,F_c,theta,T_new)
    real(real64),dimension(:,:,:),intent(inout) :: T
    integer(int64), intent(in) :: h,i,j
    real(real64), intent(in) :: d_theta,d_phi,s_notrns,s_dfsn,theta
    real(real64),dimension(:,:), intent(in) :: F_c
    real(real64),dimension(:,:,:),intent(out) :: T_new

    if(i==1) then
        ! T(h,i-1,j) = T(h,i,calculate_new_lon(j))
        if(j==1) then
             ! j-1 = N_LONS
            T_new(h,i,j) =  s_notrns*(-F_c(i,j)) + &
                        s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2 - &
                        tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
                         T(h,i,j)
        elseif(j==N_LONS) then
            ! j+1 = 1
            T_new(h,i,j) =  s_notrns*(-F_c(i,j)) + &
                        s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2 - &
                        tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                         T(h,i,j)
        else
            T_new(h,i,j) =  s_notrns*(-F_c(i,j)) + &
                        s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2 - &
                        tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                         T(h,i,j)
        end if

    elseif(i==N_LATS) then
        !T(h,i+1,j) = T(h,i,calculate_new_lon(j))
        if(j==1) then
            ! j-1 = N_LONS
            T_new(h,i,j) =  s_notrns*(-F_c(i,j)) + &
                        s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2 - &
                        tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
                         T(h,i,j)
        elseif(j==N_LONS) then
            ! j+1 = 1
            T_new(h,i,j) =  s_notrns*(-F_c(i,j)) + &
                        s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2 - &
                        tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                         T(h,i,j)
        else

            T_new(h,i,j) =  s_notrns*(-F_c(i,j)) + &
                        s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2 - &
                        tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                         T(h,i,j)
        end if

    else
        if(j==1) then
            ! j-1 = N_LONS
            T_new(h,i,j) =  s_notrns*(-F_c(i,j)) + &
                        s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2 - &
                        tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
                         T(h,i,j)
        elseif(j==N_LONS) then
            ! j+1 = 1
            T_new(h,i,j) =  s_notrns*(-F_c(i,j)) + &
                        s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2 - &
                        tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                         T(h,i,j)
        else

            T_new(h,i,j) =  s_notrns*(-F_c(i,j)) + &
                        s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2 - &
                        tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)+ &
                        1/cos(theta)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                         T(h,i,j)
        end if

    end if

end subroutine calc_new_T_deep_diff

!subroutine calc_new_T_surf_ekman(T,i,j,h,d_theta,d_phi,s_notrns,s_dfsn,F_c,theta,T_new)
!
!end subroutine


end module
