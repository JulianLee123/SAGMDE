import java.util.*;
import java.lang.Math;

class CenterComparator implements Comparator<Center>{
	//Sort centers in ascending order
    public int compare(Center s1, Center s2) {
      float d1 = 0, d2 = 0;
      for(int i = 0; i < s1.cLoc.length; i++){
        d1 += Math.pow(Math.abs(s1.cLoc[i]-s1.pLoc[i]),2);
      }
      for(int i = 0; i < s2.cLoc.length; i++){
        d2 += Math.pow(Math.abs(s2.cLoc[i]-s2.pLoc[i]),2);
      }
      if(d2 - d1 > 0)
        return -1;
      return 1;
    }
}