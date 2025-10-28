/**
Here I copy the slice output function to Vslice function -> when output the vslice, we can choose only the 80m below
*/
double THE_ang = M_PI * 55. / 180.;

struct sOutputSlice
{
    scalar *list;
    FILE *fp;
    int n;
    bool linear;
    coord plane;
};

trace void output_slice(struct sOutputSlice p)
{
    if (!p.list)
        p.list = all;
    if (p.n == 0)
        p.n = N;
    if (!p.fp)
        p.fp = stdout;
    if (!p.plane.x)
        p.plane.x = 1.;
    if (!p.plane.y)
        p.plane.y = 0.;
    if (!p.plane.z)
        p.plane.z = 1.;
    p.n++;

    int len = list_len(p.list);
    double **field = (double **)matrix_new(p.n, p.n, len * sizeof(double));
    double Delta = 0.999999 * L0 / (p.n - 1);

    // the loop of resolution grid
    // find the cartesion grid coordinate
    for (int i = 0; i < p.n; i++)
    {
        double varCoord1 = Delta * i; // some clever way of implementing general variation of coordinates instead of mapping them directly to x, y or z
        bool varX = !(p.plane.x < 1.);
        double x = (!varX ? p.plane.x * L0 : varCoord1) + X0;

        for (int j = 0; j < p.n; j++)
        {
            double varCoord2 = Delta * j;
            double y = (varX ? (p.plane.y < 1. ? p.plane.y * L0 : varCoord2) : varCoord1) + Y0;
            double z = (p.plane.z < 1. ? p.plane.z * L0 : varCoord2) + Z0;
            if (p.linear)
            {
                int k = 0;
                for (scalar s in p.list)
                    field[i][len * j + k++] = interpolate(s, x, y, z);
            }
            else
            {
                Point point = locate(x, y, z);
                int k = 0;
                for (scalar s in p.list)
                    field[i][len * j + k++] = point.level >= 0 ? s[] : nodata;
            }
        }
    }

    if (pid() == 0)
    { // master
        @ if _MPI
            MPI_Reduce(MPI_IN_PLACE, field[0], len * p.n * p.n, MPI_DOUBLE, MPI_MIN, 0,
                       MPI_COMM_WORLD);
        @endif

            // fprintf (p.fp, "x=%g\ty=%g\tz=%g\tn=%d\tlen=%d\n", p.plane.x*L0, p.plane.y*L0, p.plane.z*L0, p.n, len);
            int k = 0;
        for (scalar s in p.list)
        {
            // fprintf (p.fp, "%s\n", s.name);
            for (int i = 0; i < p.n; i++)
            {
                for (int j = 0; j < p.n; j++)
                {
                    //   fprintf(p.fp, "%g\t", (float)field[i][len * j + k]);
                    fwrite(&field[i][len * j + k], sizeof(double), 1, p.fp);
                }
                // fputc('\n', p.fp);
            }
            k++;
        }
        fflush(p.fp);
    }
    @ if _MPI 
    else // slave
        MPI_Reduce(field[0], NULL, len * p.n * p.n, MPI_DOUBLE, MPI_MIN, 0,
                   MPI_COMM_WORLD);
    @endif
    matrix_free(field);
}

trace void output_Vslice(struct sOutputSlice p)
{
    if (!p.list)
        p.list = all;
    if (p.n == 0)
        p.n = N;
    if (!p.fp)
        p.fp = stdout;
    if (!p.plane.x)
        p.plane.x = 1.;
    if (!p.plane.y)
        p.plane.y = 1.;
    if (!p.plane.z)
        p.plane.z = 0.;
    p.n++;

    int len = list_len(p.list);
    double **field = (double **)matrix_new(p.n, p.n, len * sizeof(double));
    double Delta = 0.999999 * L0 / (p.n - 1);

    // the loop of resolution grid
    // find the cartesion grid coordinate
    for (int i = 0; i < p.n; i++)
    {
        double varCoord1 = Delta * i; // some clever way of implementing general variation of coordinates instead of mapping them directly to x, y or z
        bool varX = !(p.plane.x < 1.);
        double x = (!varX ? p.plane.x * L0 : varCoord1) + X0;

        for (int j = 0; j < p.n; j++)
        {
            double varCoord2 = Delta * j;
            double y = (varX ? (p.plane.y < 1. ? p.plane.y * L0 : varCoord2) : varCoord1) + Y0;
            double z = (p.plane.z < 1. ? p.plane.z * L0 : varCoord2) + Z0;
            if (p.linear)
            {
                int k = 0;
                for (scalar s in p.list)
                    field[i][len * j + k++] = interpolate(s, x, y, z);
            }
            else
            {
                Point point = locate(x, y, z);
                int k = 0;
                for (scalar s in p.list)
                    field[i][len * j + k++] = point.level >= 0 ? s[] : nodata;
            }
        }
    }

    if (pid() == 0)
    { // master
        @ if _MPI
            MPI_Reduce(MPI_IN_PLACE, field[0], len * p.n * p.n, MPI_DOUBLE, MPI_MIN, 0,
                       MPI_COMM_WORLD);
        @endif

            // fprintf (p.fp, "x=%g\ty=%g\tz=%g\tn=%d\tlen=%d\n", p.plane.x*L0, p.plane.y*L0, p.plane.z*L0, p.n, len);
            int k = 0;
        for (scalar s in p.list)
        {
            // fprintf (p.fp, "%s\n", s.name);
            for (int i = 0; i < p.n; i++)
            {
                for (int j = 0; j < 100; j++)
                {
                    // fprintf(p.fp, "%g\t", (float)field[i][len * j + k]);
                    fwrite(&field[i][len * j + k], sizeof(double), 1, p.fp);
                }
                // fputc('\n', p.fp);
            }
            k++;
        }
        fflush(p.fp);
    }
    @ if _MPI 
    else // slave
        MPI_Reduce(field[0], NULL, len * p.n * p.n, MPI_DOUBLE, MPI_MIN, 0,
                   MPI_COMM_WORLD);
    @endif
    matrix_free(field);
}

trace 
void output_Vyslice(struct sOutputSlice p)
{
    if (!p.list)
        p.list = all;
    if (p.n == 0)
        p.n = N;
    if (!p.fp)
        p.fp = stdout;
    if (!p.plane.x)
        p.plane.x = 1.;
    if (!p.plane.y)
        p.plane.y = 1.;
    if (!p.plane.z)
        p.plane.z = 0.;
    p.n++;

    int len = list_len(p.list);
    double **field = (double **)matrix_new(p.n, p.n, len * sizeof(double));
    double Delta = 0.999999 * L0 / (p.n - 1);

    // the loop of resolution grid
    // find the cartesion grid coordinate
    for (int i = 0; i < p.n; i++)
    {
        double varCoord1 = Delta * i; // some clever way of implementing general variation of coordinates instead of mapping them directly to x, y or z
        bool varX = !(p.plane.x < 1.);
        double x = (!varX ? p.plane.x * L0 : varCoord1) + X0;

        for (int j = 0; j < p.n; j++)
        {
            double varCoord2 = Delta * j;
            double y = (varX ? (p.plane.y < 1. ? p.plane.y * L0 : varCoord2) : varCoord1) + Y0;
            double z = (p.plane.z < 1. ? p.plane.z * L0 : varCoord2) + Z0;
            if (p.linear)
            {
                int k = 0;
                for (scalar s in p.list)
                    field[i][len * j + k++] = interpolate(s, x, y, z);
            }
            else
            {
                Point point = locate(x, y, z);
                int k = 0;
                for (scalar s in p.list)
                    field[i][len * j + k++] = point.level >= 0 ? s[] : nodata;
            }
        }
    }

    if (pid() == 0)
    { // master
        @ if _MPI
            MPI_Reduce(MPI_IN_PLACE, field[0], len * p.n * p.n, MPI_DOUBLE, MPI_MIN, 0,
                       MPI_COMM_WORLD);
        @endif

            // fprintf (p.fp, "x=%g\ty=%g\tz=%g\tn=%d\tlen=%d\n", p.plane.x*L0, p.plane.y*L0, p.plane.z*L0, p.n, len);
            int k = 0;
        for (scalar s in p.list)
        {
            // fprintf (p.fp, "%s\n", s.name);
            for (int i = 500; i < 520; i++)
            {
                for (int j = 0; j < p.n; j++)
                {
                    // fprintf(p.fp, "%g\t", (float)field[i][len * j + k]);
                    fwrite(&field[i][len * j + k], sizeof(double), 1, p.fp);
                }
                // fputc('\n', p.fp);
            }
            k++;
        }
        fflush(p.fp);
    }
    @ if _MPI 
    else // slave
        MPI_Reduce(field[0], NULL, len * p.n * p.n, MPI_DOUBLE, MPI_MIN, 0,
                   MPI_COMM_WORLD);
    @endif
    matrix_free(field);
}


trace void output_Vslice_TV(struct sOutputSlice p)
{
    if (!p.list)
        p.list = all;
    if (p.n == 0)
        p.n = N;
    if (!p.fp)
        p.fp = stdout;
    if (!p.plane.x)
        p.plane.x = 1.;
    if (!p.plane.y)
        p.plane.y = 1.;
    if (!p.plane.z)
        p.plane.z = 0.;
    p.n++;

    int len = list_len(p.list);
    double **field = (double **)matrix_new(p.n, p.n, len * sizeof(double));
    double Delta = 0.999999 * L0 / (p.n - 1);

    // the loop of resolution grid
    // find the cartesion grid coordinate
    for (int i = 0; i < p.n; i++)
    {
        double varCoord1 = Delta * i; // some clever way of implementing general variation of coordinates instead of mapping them directly to x, y or z
        bool varX = !(p.plane.x < 1.);
        double x = (!varX ? p.plane.x * L0 : varCoord1) + X0;

        for (int j = 0; j < p.n; j++)
        {
            double varCoord2 = Delta * j;
            double y = (varX ? (p.plane.y < 1. ? p.plane.y * L0 : varCoord2) : varCoord1) + Y0;
            double z = (p.plane.z < 1. ? p.plane.z * L0 : varCoord2) + Z0;
            if (p.linear)
            {
                int k = 0;
                for (scalar s in p.list)
                    field[i][len * j + k++] = interpolate(s, x, y, z);
            }
            else
            {
                Point point = locate(x, y, z);
                int k = 0;
                for (scalar s in p.list)
                    field[i][len * j + k++] = point.level >= 0 ? s[] : nodata;
            }
        }
    }

    if (pid() == 0)
    { // master
        @ if _MPI
            MPI_Reduce(MPI_IN_PLACE, field[0], len * p.n * p.n, MPI_DOUBLE, MPI_MIN, 0,
                       MPI_COMM_WORLD);
        @endif

            // fprintf (p.fp, "x=%g\ty=%g\tz=%g\tn=%d\tlen=%d\n", p.plane.x*L0, p.plane.y*L0, p.plane.z*L0, p.n, len);
            int k = 0;
        for (scalar s in p.list)
        {
            // fprintf (p.fp, "%s\n", s.name);
            for (int i = 0; i < p.n; i++)
            {
                for (int j = 0; j < 5; j++)
                {
                    // fprintf(p.fp, "%g\t", (float)field[i][len * j + k]);
                    fwrite(&field[i][len * j + k], sizeof(double), 1, p.fp);
                }
                // fputc('\n', p.fp);
            }
            k++;
        }
        fflush(p.fp);
    }
    @ if _MPI 
    else // slave
        MPI_Reduce(field[0], NULL, len * p.n * p.n, MPI_DOUBLE, MPI_MIN, 0,
                   MPI_COMM_WORLD);
    @endif
    matrix_free(field);
}

trace void output_W12(struct sOutputSlice p)
{
    if (!p.list)
        p.list = all;
    if (p.n == 0)
        p.n = N;
    if (!p.fp)
        p.fp = stdout;
    if (!p.plane.x)
        p.plane.x = 1.;
    if (!p.plane.y)
        p.plane.y = 1.;
    if (!p.plane.z)
        p.plane.z = 0.;
    // p.n++;

    int len = list_len(p.list);
    double **field = (double **)matrix_new(p.n, p.n, len * sizeof(double));
    double Delta = 0.999999 * L0 / 2048.;

    // the W12 output (11 * 11): i = y (0..10),  j = z (990..1100), x = 528
    for (int i = 0; i < p.n; i++)
    {
        double varCoord1 = Delta * i; // some clever way of implementing general variation of coordinates instead of mapping them directly to x, y or z
        bool varX = !(p.plane.x < 1.);
        double x = (!varX ? p.plane.x * L0 : varCoord1) + X0;

        for (int j = 0; j < p.n; j++)
        {
            double varCoord2 = Delta * (j + 990.);
            double y = (varX ? (p.plane.y < 1. ? p.plane.y * L0 : varCoord2) : varCoord1) + Y0;
            double z = (p.plane.z < 1. ? p.plane.z * L0 : varCoord2) + Z0;
            if (p.linear)
            {
                int k = 0;
                for (scalar s in p.list)
                    field[i][len * j + k++] = interpolate(s, x, y, z);
            }
            else
            {
                Point point = locate(x, y, z);
                int k = 0;
                for (scalar s in p.list)
                    field[i][len * j + k++] = point.level >= 0 ? s[] : nodata;
            }
        }
    }

    if (pid() == 0)
    { // master
        @ if _MPI
            MPI_Reduce(MPI_IN_PLACE, field[0], len * p.n * p.n, MPI_DOUBLE, MPI_MIN, 0,
                       MPI_COMM_WORLD);
        @endif

            // fprintf (p.fp, "x=%g\ty=%g\tz=%g\tn=%d\tlen=%d\n", p.plane.x*L0, p.plane.y*L0, p.plane.z*L0, p.n, len);
            int k = 0;
        for (scalar s in p.list)
        {
            // fprintf (p.fp, "%s\n", s.name);
            for (int i = 0; i < p.n; i++)
            {
                for (int j = 0; j < p.n; j++)
                {
                    // fprintf(p.fp, "%g\t", (float)field[i][len * j + k]);
                    fwrite(&field[i][len * j + k], sizeof(double), 1, p.fp);
                }
                // fputc('\n', p.fp);
            }
            k++;
        }
        fflush(p.fp);
    }
    @ if _MPI 
    else // slave
        MPI_Reduce(field[0], NULL, len * p.n * p.n, MPI_DOUBLE, MPI_MIN, 0,
                   MPI_COMM_WORLD);
    @endif
    matrix_free(field);
}


trace void b12output_slice(struct sOutputSlice p)
{
    if (!p.list)
        p.list = all;
    if (p.n == 0)
        p.n = N;
    if (!p.fp)
        p.fp = stdout;
    if (!p.plane.x)
        p.plane.x = 1.;
    if (!p.plane.y)
        p.plane.y = 1.;
    if (!p.plane.z)
        p.plane.z = 1.;
    p.n++;

    int len = list_len(p.list);
    double **field = (double **)matrix_new(p.n, p.n, len * sizeof(double));
    int H = p.plane.y;

    // the loop of resolution grid
    // find the cartesion grid coordinate
    for (int i = 0; i < p.n; i++)
    {
        double NS_hist = i - p.n / 2.;
        double xf0 = 500. + NS_hist * cos(THE_ang);
        double zf0 = L0 / 2. + NS_hist * sin(THE_ang);
        for (int j = 0; j < p.n; j++)
        {
            double WE_dist = j - p.n / 2.;
            double xx = xf0 - WE_dist * sin(THE_ang);
            double zz = zf0 + WE_dist * cos(THE_ang);
            if (p.linear)
            {
                int k = 0;
                for (scalar s in p.list)
                    field[i][len * j + k++] = interpolate(s, xx, H, zz);
            }
            else
            {
                Point point = locate(xx, H, zz);
                int k = 0;
                for (scalar s in p.list)
                    field[i][len * j + k++] = point.level >= 0 ? s[] : nodata;
            }
        }
    }

    if (pid() == 0)
    { // master
        @ if _MPI
            MPI_Reduce(MPI_IN_PLACE, field[0], len * p.n * p.n, MPI_DOUBLE, MPI_MIN, 0,
                       MPI_COMM_WORLD);
        @endif

            // fprintf (p.fp, "x=%g\ty=%g\tz=%g\tn=%d\tlen=%d\n", p.plane.x*L0, p.plane.y*L0, p.plane.z*L0, p.n, len);
            int k = 0;
        for (scalar s in p.list)
        {
            // fprintf (p.fp, "%s\n", s.name);
            for (int i = 0; i < p.n; i++)
            {
                for (int j = 0; j < p.n; j++)
                {
                    fprintf(p.fp, "%g\t", (float)field[i][len * j + k]);
                }
                fputc('\n', p.fp);
            }
            k++;
        }
        fflush(p.fp);
    }
    @ if _MPI 
    else // slave
        MPI_Reduce(field[0], NULL, len * p.n * p.n, MPI_DOUBLE, MPI_MIN, 0,
                   MPI_COMM_WORLD);
    @endif
    matrix_free(field);
}