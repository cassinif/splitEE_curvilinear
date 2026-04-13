function T = act_3(T,V,iV,fd1d2d3)
  T = mump(fd1d2d3.*mump(T,iV,3),V,3);
end
