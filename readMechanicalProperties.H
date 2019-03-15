    Info<< "Reading mechanical properties\n" << endl;

    IOdictionary mechanicalProperties
    (
        IOobject
        (
            "mechanicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const dictionary& rhoDict(mechanicalProperties.subDict("rho"));
    word rhoType(rhoDict.lookup("type"));

    autoPtr<volScalarField> rhoPtr;

    IOobject rhoIO
    (
        "rho",
        runTime.timeName(0),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (rhoType == "uniform")
    {
        scalar rhoValue(readScalar(rhoDict.lookup("value")));

        rhoPtr.reset
        (
            new volScalarField
            (
                rhoIO,
                mesh,
                dimensionedScalar
                (
                    "rho",
                    dimMass/dimVolume,
                    rhoValue
                )
            )
        );
    }
    else if (rhoType == "field")
    {
        rhoIO.readOpt() = IOobject::MUST_READ;

        rhoPtr.reset
        (
            new volScalarField
            (
                rhoIO,
                mesh
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are uniform or field for rho"
            << abort(FatalError);
    }

    volScalarField& rho = rhoPtr();

// START ALID
    Info<< "Setting up ALID\n" << endl;
    const dictionary& alidDict(mechanicalProperties.subDict("ALID"));
    word alidType(alidDict.lookup("type"));

    autoPtr<volScalarField> alidPtr;

    IOobject alidIO
    (
        "ALID",
        runTime.timeName(0),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (alidType == "field")
    {
        alidIO.readOpt() = IOobject::MUST_READ;

        alidPtr.reset
        (
            new volScalarField
            (
                alidIO,
                mesh
            )
        );
    }
    else if (alidType == "none" || alidType =="uniform")
    {
        scalar alidValue(0.0);
        
	    alidPtr.reset
        (
            new volScalarField
            (
                alidIO,
                mesh,
                dimensionedScalar
                (
                    "ALID",
                    dimless/dimTime,
                    alidValue
                )
            )
        );

    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are uniform, field or none for ALID"
            << abort(FatalError);
    }

    volScalarField& alid = alidPtr();

    if (alidType == "uniform")
    {
        scalar Rstart(readScalar(alidDict.lookup("Rstart")));
        scalar Rend(readScalar(alidDict.lookup("Rend")));
        scalar Cmax(readScalar(alidDict.lookup("Cmax")));

	    scalar dR = Rend - Rstart; 
	    forAll(mesh.C(), celli)
	    {
	        scalar r_ = 0.0;
	        for (int j=0; j<3; j++)
	        {
	            r_ += mesh.C()[celli][j] * mesh.C()[celli][j];
 	        }
	        r_ = Foam::sqrt(r_);

	        if (r_>=Rstart && r_<=Rend)
	        {
		        alid[celli] = Foam::pow( (r_-Rstart) / dR, 3.0 );
	        }
	        if (r_>Rend)
	        {
		        alid[celli] = 1.0;
  	        }
	        alid[celli] = Cmax * alid[celli];
	    }
    }

// END ALID

    const dictionary& EDict(mechanicalProperties.subDict("E"));
    word EType(EDict.lookup("type"));

    autoPtr<volScalarField> EPtr;

    IOobject EHeader
    (
        "E",
        runTime.timeName(0),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (EType == "uniform")
    {
        scalar rhoEValue(readScalar(EDict.lookup("value")));

        EPtr.reset
        (
            new volScalarField
            (
                EHeader,
                mesh,
                dimensionedScalar
                (
                    "Erho",
                    dimMass/dimLength/sqr(dimTime),
                    rhoEValue
                )
            )
        );
    }
    else if (EType == "field")
    {
        EHeader.readOpt() = IOobject::MUST_READ;

        EPtr.reset
        (
            new volScalarField
            (
                EHeader,
                mesh
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are uniform or field for E"
            << abort(FatalError);
    }

    volScalarField& rhoE = EPtr();

    autoPtr<volScalarField> nuPtr;

    IOobject nuIO
    (
        "nu",
        runTime.timeName(0),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    const dictionary& nuDict(mechanicalProperties.subDict("nu"));
    word nuType(nuDict.lookup("type"));

    if (nuType == "uniform")
    {
        scalar nuValue(readScalar(nuDict.lookup("value")));
        nuPtr.reset
        (
            new volScalarField
            (
                nuIO,
                mesh,
                dimensionedScalar
                (
                    "nu",
                    dimless,
                    nuValue
                )
            )
        );
    }
    else if (nuType == "field")
    {
        nuIO.readOpt() = IOobject::MUST_READ;
        nuPtr.reset
        (
            new volScalarField
            (
                nuIO,
                mesh
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are uniform or field for nu"
            << abort(FatalError);
    }

    volScalarField& nu = nuPtr();

    Info<< "Normalising E : E/rho\n" << endl;
    volScalarField E(rhoE/rho);

    Info<< "Calculating Lame's coefficients\n" << endl;

    volScalarField mu(E/(2.0*(1.0 + nu)));
    volScalarField lambda(nu*E/((1.0 + nu)*(1.0 - 2.0*nu)));
    volScalarField threeK(E/(1.0 - 2.0*nu));

    if (mechanicalProperties.get<bool>("planeStress"))
    {
        Info<< "Plane Stress\n" << endl;

        lambda = nu*E/((1.0 + nu)*(1.0 - nu));
        threeK = E/(1.0 - nu);
    }
    else
    {
        Info<< "Plane Strain\n" << endl;
    }
