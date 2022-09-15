!********************************************************************************************************
! Module contains subroutines to calculate the new temperature at new time step using Forward Euler method
!********************************************************************************************************
module calc_new_T_fe
use Constants
use DEGREE_TO_RADIAN
use dA_da
use div_m
use, intrinsic :: iso_fortran_env
implicit none
public:: calc_new_T_surf_diff, calc_new_T_deep_diff,calc_new_T_surf_notrns, &
        calc_new_T_deep_notrns,calc_new_T_surf_ekman,calc_new_T_deep_ekman
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

!***************************************************
!Function to calculate new coordinates at boundaries
!***************************************************
function calculate_new_lon(lon) result(new_lon)
! At maximum or minimum latitude, the adjacent latitude point above or bellow, respectively, is located at the same latitude,
! but at the new longitudinal point (from lon to new_lon)
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

        !print*, 'theta=',theta
        ! T(h,i-1,j) = T(h,i,calculate_new_lon(j))
        if(j==1) then
             ! j-1 = N_LONS
            T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
            s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2  - &
            tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)  + &
            (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
            T(h,i,j)

        elseif(j==N_LONS) then
            ! j+1 = 1
            T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
            s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2  - &
            tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)  + &
            (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
            T(h,i,j)

        else
            T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
            s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2 - &
            tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta) + &
            (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
            T(h,i,j)


        end if

    elseif(i==N_LATS) then
        !T(h,i+1,j) = T(h,i,calculate_new_lon(j))

        if(j==1) then
            ! j-1 = N_LONS
            T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
            s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
            tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)   + &
            (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
            T(h,i,j)
        elseif(j==N_LONS) then
            ! j+1 = 1
            T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
            s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
            tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)  + &
            (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
            T(h,i,j)
        else

            T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
            s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
            tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)  + &
            (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
            T(h,i,j)
        end if

    else

        if(j==1) then
            ! j-1 = N_LONS
            T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
            s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
            tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
            (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
            T(h,i,j)
        elseif(j==N_LONS) then
            ! j+1 = 1
            T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
            s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
            tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
            (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
            T(h,i,j)
        else

            T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
            s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
            tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
            (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
            T(h,i,j)
        end if

    endif


end subroutine calc_new_T_surf_diff

subroutine calc_new_T_deep_diff(T,i,j,h,d_theta,d_phi,s_notrns,s_dfsn,F_c,theta,T_new)
    real(real64),dimension(:,:,:),intent(in) :: T
    integer(int64), intent(in) :: h,i,j
    real(real64), intent(in) :: d_theta,d_phi,s_notrns,s_dfsn,theta
    real(real64),dimension(:,:), intent(in) :: F_c
    real(real64),dimension(:,:,:),intent(out) :: T_new

    if(i==1) then
        ! T(h,i-1,j) = T(h,i,calculate_new_lon(j))

        if(j==1) then
             ! j-1 = N_LONS
            T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
            s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2  - &
            tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)  + &
            (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
            T(h,i,j)
        elseif(j==N_LONS) then
            ! j+1 = 1
            T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
            s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2  - &
            tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)  + &
            (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
            T(h,i,j)
        else
            T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
            s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2  - &
            tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)  + &
            (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
            T(h,i,j)

        end if

    elseif(i==N_LATS) then

        !T(h,i+1,j) = T(h,i,calculate_new_lon(j))
        if(j==1) then
            ! j-1 = N_LONS
            T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
            s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
            tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)  + &
            (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
            T(h,i,j)
        elseif(j==N_LONS) then
            ! j+1 = 1
            T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
            s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
            tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)  + &
            (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
            T(h,i,j)
        else

            T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
            s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
            tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)  + &
            (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
            T(h,i,j)
        end if

    else

        if(j==1) then
            ! j-1 = N_LONS
            T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
            s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2 - &
            tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
            (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
            T(h,i,j)
        elseif(j==N_LONS) then
            ! j+1 = 1
            T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
            s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
            tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
            (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
            T(h,i,j)
        else
            T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
            s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
            tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
            (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
            T(h,i,j)
        end if

    end if

end subroutine calc_new_T_deep_diff

!*************************************************************************************************************
! Subroutines to calculate new temperature using the forward Euler method (with diffusion and Ekman transport)
!*************************************************************************************************************
subroutine calc_new_T_surf_ekman(T,M,i,j,h,d_theta,d_phi,s_notrns,s_dfsn,s_ekmn,F_c,F_a,theta,T_new)
    real(real64),dimension(:,:,:),intent(in) :: T,M
    real(real64), intent(in) :: d_theta,d_phi,s_notrns,s_dfsn,s_ekmn,theta
    integer(int64),intent(in) :: i,j,h
    real(real64),dimension(:,:),intent(in) :: F_c,F_a
    real(real64) :: div_m,dM_theta_dtheta,dT_surf_dtheta,dT_surf_dphi,dM_phi_dphi
    real(real64),dimension(:,:,:), intent(out) :: T_new

    !Indexes for Theta and PHI: Theta = 1 , PHI =2 .

    call calculate_div_M(i,j,theta,M(2,:,:),M(1,:,:),div_m)
    call calculate_dA_d_phi(T(1,:,:),i,j,dT_surf_dphi)
    call calculate_dA_d_phi(M(2,:,:),i,j,dM_phi_dphi)
    call calculate_dA_d_theta(T(1,:,:),i,j,dT_surf_dtheta)
    call calculate_dA_d_theta(M(1,:,:),i,j,dM_theta_dtheta)

    if(div_m < 0) then

        if(i==1) then
            ! T(h,i-1,j) = T(h,i,calculate_new_lon(j))
            if(j==1) then
                 ! j-1 = N_LONS
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            elseif(j==N_LONS) then
                ! j+1 = 1
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            else
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2 - &
                tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta) + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)



            end if

        elseif(i==N_LATS) then
            !T(h,i+1,j) = T(h,i,calculate_new_lon(j))
            if(j==1) then
                ! j-1 = N_LONS
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)   + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            elseif(j==N_LONS) then
                ! j+1 = 1
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            else
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)


            end if

        else

            if(j==1) then
                ! j-1 = N_LONS
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)


            elseif(j==N_LONS) then
                ! j+1 = 1
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)


            else
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            end if

        end if

    elseif(div_m > 0) then

        if(i==1) then
            ! T(h,i-1,j) = T(h,i,calculate_new_lon(j))
            if(j==1) then
                 ! j-1 = N_LONS
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) - &
                 s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            elseif(j==N_LONS) then
                ! j+1 = 1
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            else
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2 - &
                tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta) + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            end if

        elseif(i==N_LATS) then
            !T(h,i+1,j) = T(h,i,calculate_new_lon(j))
            if(j==1) then
                ! j-1 = N_LONS
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)   + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            elseif(j==N_LONS) then
                ! j+1 = 1
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            else
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            end if

        else

            if(j==1) then
                ! j-1 = N_LONS
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            elseif(j==N_LONS) then
                ! j+1 = 1
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            else
                T_new(h,i,j) = s_notrns*(F_a(i,j)+F_c(i,j)-SIGMA*EMISSIVITY*(T(h,i,j))**4)+ &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) - &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_surf_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_surf_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            end if

        endif


    end if

end subroutine

subroutine calc_new_T_deep_ekman(T,M,i,j,h,d_theta,d_phi,s_notrns,s_dfsn,s_ekmn,F_c,theta,T_new)

    real(real64),dimension(:,:,:),intent(in) :: T,M
    real(real64), intent(in) :: d_theta,d_phi,s_notrns,s_dfsn,s_ekmn,theta
    integer(int64),intent(in) :: i,j,h
    real(real64),dimension(:,:),intent(in) :: F_c
    real(real64) :: div_m,dM_theta_dtheta,dT_deep_dtheta,dT_deep_dphi,dM_phi_dphi
    real(real64),dimension(:,:,:), intent(out) :: T_new

    call calculate_div_M(i,j,theta,M(2,:,:),M(1,:,:),div_m) ! remember div_m is multiplied by 1/R_planet here
    call calculate_dA_d_phi(T(2,:,:),i,j,dT_deep_dphi)
    call calculate_dA_d_phi(M(2,:,:),i,j,dM_phi_dphi)
    call calculate_dA_d_theta(T(2,:,:),i,j,dT_deep_dtheta)
    call calculate_dA_d_theta(M(1,:,:),i,j,dM_theta_dtheta)

    if(div_m < 0) then
        if(i==1) then
            ! T(h,i-1,j) = T(h,i,calculate_new_lon(j))

            if(j==1) then
                 ! j-1 = N_LONS
                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            elseif(j==N_LONS) then
                ! j+1 = 1
                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            else
                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            end if

        elseif(i==N_LATS) then

            !T(h,i+1,j) = T(h,i,calculate_new_lon(j))
            if(j==1) then
                ! j-1 = N_LONS
                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            elseif(j==N_LONS) then
                ! j+1 = 1
                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            else

                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            end if

        else

            if(j==1) then
                ! j-1 = N_LONS
                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2 - &
                tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            elseif(j==N_LONS) then
                ! j+1 = 1
                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            else
                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(1,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            end if

        end if

    elseif(div_m > 0) then
        if(i==1) then
            ! T(h,i-1,j) = T(h,i,calculate_new_lon(j))

            if(j==1) then
                 ! j-1 = N_LONS
                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            elseif(j==N_LONS) then
                ! j+1 = 1
                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)
            else
                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i,calculate_new_lon(j)))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i,calculate_new_lon(j)))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            end if

        elseif(i==N_LATS) then

            !T(h,i+1,j) = T(h,i,calculate_new_lon(j))
            if(j==1) then
                ! j-1 = N_LONS
                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)
            elseif(j==N_LONS) then
                ! j+1 = 1
                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)
            else

                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i,calculate_new_lon(j))-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i,calculate_new_lon(j))-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            end if

        else

            if(j==1) then
                ! j-1 = N_LONS
                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2 - &
                tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,N_LONS))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            elseif(j==N_LONS) then
                ! j+1 = 1
                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            else
                T_new(h,i,j) = s_notrns*(-F_c(i,j)) + &
                s_dfsn*( (T(h,i+1,j)-2*T(h,i,j)+T(h,i-1,j))/(d_theta)**2  - &
                tan(theta)*(T(h,i+1,j)-T(h,i-1,j))/(2*d_theta)  + &
                (1/(cos(theta))**2)*(T(h,i,j+1)-2*T(h,i,j)+T(h,i,j-1))/(d_phi**2)) + &
                s_ekmn*(T(h,i,j)*dM_theta_dtheta + M(1,i,j)*dT_deep_dtheta-T(h,i,j)*M(1,i,j)*tan(theta)+ &
                (1/cos(theta))*(M(2,i,j)*dT_deep_dphi+T(h,i,j)*dM_phi_dphi) - &
                T(2,i,j)*(dM_theta_dtheta-M(1,i,j)*tan(theta)+1/cos(theta)*dM_phi_dphi) ) + &
                T(h,i,j)

            end if

        end if

    end if


end subroutine


















end module
